#ifdef GEOMETRY_DEBUG

#ifdef G4VIS_USE_OPENGLX
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLStoredX.hh"
#endif

#ifdef G4VIS_USE_OPENGLXM
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLStoredXm.hh"
#endif

#include "G4ASCIITree.hh"
#include "PadVisManager.hh"

PadVisManager::PadVisManager() {;}

PadVisManager::~PadVisManager() {;}

void PadVisManager::RegisterGraphicsSystems()
{
#ifdef G4VIS_USE_OPENGLX
  RegisterGraphicsSystem(new G4OpenGLImmediateX);
  RegisterGraphicsSystem(new G4OpenGLStoredX);
#endif

#ifdef G4VIS_USE_OPENGLXM
  RegisterGraphicsSystem(new G4OpenGLImmediateXm);
  RegisterGraphicsSystem(new G4OpenGLStoredXm);
#endif

  RegisterGraphicsSystem(new G4ASCIITree);

    G4cout <<
      "\nYou have successfully chosen to use the following graphics systems."
	   << G4endl;
    //PrintAvailableGraphicsSystems();
}

#endif
