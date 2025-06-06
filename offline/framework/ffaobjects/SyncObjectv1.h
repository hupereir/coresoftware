// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_SYNCOBJECTV1_H
#define FFAOBJECTS_SYNCOBJECTV1_H

#include "SyncObject.h"

#include <iostream>

class PHObject;

class SyncObjectv1 : public SyncObject
{
 public:
  /// ctor
  SyncObjectv1() = default;
  explicit SyncObjectv1(const SyncObject& source);

  PHObject* CloneMe() const override { return new SyncObjectv1(*this); }
  /// dtor
  ~SyncObjectv1() override = default;

  ///  Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& out = std::cout) const override;

  /// isValid returns non zero if object contains valid data
  int isValid() const override;

  /// set Event Counter
  // cppcheck-suppress virtualCallInConstructor
  void EventCounter(const int ival) override { eventcounter = ival; }

  /// set Event Number
  // cppcheck-suppress virtualCallInConstructor
  void EventNumber(const int ival) override { eventnumber = ival; }

  /// set Run Number
  // cppcheck-suppress virtualCallInConstructor
  void RunNumber(const int ival) override { runnumber = ival; }

  /// set Segment Number
  // cppcheck-suppress virtualCallInConstructor
  void SegmentNumber(const int ival) override { segmentnumber = ival; }

 protected:
  /// get Event Counter
  int EventCounter() const override { return eventcounter; }
  /// get Event Number
  int EventNumber() const override { return eventnumber; }
  /// get Run Number
  int RunNumber() const override { return runnumber; }
  /// get Segment Number
  int SegmentNumber() const override { return segmentnumber; }

 private:
  int eventcounter = 0;         // running counter
  int eventnumber = 0;          // Event number
  int runnumber = 0;            // Run number
  int segmentnumber = -999999;  // segment number

  ClassDefOverride(SyncObjectv1, 1)
};

#endif
