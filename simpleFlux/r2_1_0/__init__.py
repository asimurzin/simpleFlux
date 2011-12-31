#!/usr/bin/env python

#---------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#---------------------------------------------------------------------------
from Foam import man, ref


#---------------------------------------------------------------------------
def createFields( runTime, mesh ):

  ref.ext_Info() << "Reading field p\n" << ref.nl
  
  p = man.volScalarField( man.IOobject( ref.word( "p" ), 
                                        ref.fileName( runTime.timeName() ), 
                                        mesh, 
                                        ref.IOobject.MUST_READ, 
                                        ref.IOobject.AUTO_WRITE ), mesh )

  ref.ext_Info() << "Reading field U\n" << ref.nl
  U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.MUST_READ,
                                            ref.IOobject.AUTO_WRITE ),
                            mesh );
  
  phi = man.createPhi( runTime, mesh, U );
  
  pRefCell = 0
  pRefValue = 0.0

  pRefCell, pRefValue = ref.setRefCell( p, mesh.solutionDict().subDict( ref.word( "SIMPLE" ) ), pRefCell, pRefValue )

  laminarTransport = man.singlePhaseTransportModel( U, phi )

  turbulence = man.incompressible.RASModel.New( U, phi, laminarTransport )
  
  sources = man.IObasicSourceList( mesh )
  
  return p, U, phi, pRefCell, pRefValue, laminarTransport, turbulence, sources


#---------------------------------------------------------------------------
def fun_UEqn( U, phi, turbulence, p, sources ):
  # Solve the Momentum equation
  
  UEqn = man.fvm.div( phi, U ) + man( turbulence.divDevReff( U ), man.Deps( turbulence, U ) ) == man( sources( U ), man.Deps( U ) )

  UEqn.relax()

  sources.constrain( UEqn )

  ref.solve( UEqn == -man.fvc.grad( p ) )

  return UEqn


#---------------------------------------------------------------------------
def fun_pEqn( mesh, runTime, simple, U, phi, turbulence, p, UEqn, pRefCell, pRefValue, cumulativeContErr ):
  
  p.ext_boundaryField().updateCoeffs()

  rAU = 1.0 / UEqn().A();
  U << rAU * UEqn().H() 
  
  phi << ( ref.fvc.interpolate( U, ref.word( "interpolate(HbyA)" ) ) & mesh.Sf() )
  
  ref.adjustPhi(phi, U, p)

  # Non-orthogonal pressure corrector loop
  while simple.correctNonOrthogonal():
    pEqn = ref.fvm.laplacian( rAU, p ) == ref.fvc.div( phi )

    pEqn.setReference( pRefCell, pRefValue )

    pEqn.solve()

    if simple.finalNonOrthogonalIter():
      phi -= pEqn.flux()
      pass
    pass
  cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )

  # Explicitly relax pressure for momentum corrector
  p.relax()

  # Momentum corrector
  U -= rAU * ref.fvc.grad( p )
  U.correctBoundaryConditions()
  
  return cumulativeContErr


#---------------------------------------------------------------------------
def main_standalone( argc, argv ):
  
  args = ref.setRootCase( argc, argv )
  
  runTime = man.createTime( args )
  
  mesh = man.createMesh( runTime )
    
  p, U, phi, pRefCell, pRefValue, laminarTransport, turbulence,sources = createFields( runTime, mesh )

  cumulativeContErr = ref.initContinuityErrs()
  
  simple = man.simpleControl (mesh)

  # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
  ref.ext_Info() << "\nStarting time loop\n" << ref.nl

  while simple.loop():
    ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl

    # --- Pressure-velocity SIMPLE corrector
    UEqn = fun_UEqn( U, phi, turbulence, p, sources )
    cumulativeContErr = fun_pEqn( mesh, runTime, simple, U, phi, turbulence, p, UEqn, pRefCell, pRefValue, cumulativeContErr )

    turbulence.correct()

    runTime.write()

    ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
            << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
            << ref.nl << ref.nl
    pass

  ref.ext_Info() << "End\n" << ref.nl

  import os
  return os.EX_OK


#---------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020100" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.1.0 \n "     
    


