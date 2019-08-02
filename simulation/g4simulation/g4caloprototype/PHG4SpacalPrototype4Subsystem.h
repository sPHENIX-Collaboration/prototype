// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SPACALPROTOTYPESUBSYSTEM_H
#define G4DETECTORS_PHG4SPACALPROTOTYPESUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>  // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4SpacalPrototype4Detector;
class PHG4SteppingAction;

class PHG4SpacalPrototype4Subsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  explicit PHG4SpacalPrototype4Subsystem(const std::string& name = "CEMC");

  //! destructor
  virtual ~PHG4SpacalPrototype4Subsystem(void) {}

  //! init
  /*!
   creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
   reates the stepping action and place it on the node tree, under "ACTIONS" node
   creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
   */
  int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
   get all relevant nodes from top nodes (namely hit list)
   and pass that to the stepping action
   */
  int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector(void) const;
  PHG4SteppingAction* GetSteppingAction(void) const { return steppingAction_; }

  void
  Print(const std::string& what = "ALL") const;

 private:
  //! load the default parameter to param
  void SetDefaultParameters();

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4SpacalPrototype4Detector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;
};

#endif
