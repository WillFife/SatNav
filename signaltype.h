// signaltype.h
//
// SignalType is a unique representation of a system, frequency, and other
// attributes.
//
// Copyright (c) 2016 The GRID Software Project. All rights reserved.  Use of
// this source code is governed by the license dsprx/doc/gridLicense/GRID_LICENSE.pdf

#ifndef __SIGNAL_TYPE_GRID_H
#define __SIGNAL_TYPE_GRID_H

#include <string>
#include "typedefs.h"
#include "system.h"
#include "txid.h"

// Each SignalType is a unique combination of four components:
//   - system
//   - frequency
//   - third attribute which may be a code (e.g., C/A), and/or a channel
//     (e.g., I or Q), and/or a tracking mode (e.g., CLM for combined CL and CM
//     tracking).
//   - signal group (e.g., PRIMARY or ALT1), which is a logical group of
//     signals, typically from a specific antenna.
// The first three components comprise the GenericType.  The final component
// is known as the Group.  A SignalType object's GenericType and Group can be
// accessed separately.
//
// The GenericType enumerations are constructed as SYSTEM_FREQUENCY_ATTRIBUTE.
// Depending on the third attribute, a GenericType is either primitive or
// composite.  A primitive GenericType is associated with a single code and
// channel; a composite GenericType is associated with multiple codes or
// channels. The function isPrimitive() identifies primitive GenericTypes.

class SignalType {
 public:
  enum GenericType {
    GPS_L1_CA              = 0,  // GPS L1 civil C/A code
    GPS_L2_CM              = 1,  // GPS L2 civil M code
    GPS_L2_CL              = 2,  // GPS L2 civil L code
    GPS_L2_CLM             = 3,  // GPS L2 M+L combined code
    CDMA_UHF_PILOT_I       = 5,  // Cellular CDMA pilot code on I channel
    CDMA_UHF_PILOT_Q       = 6,  // Cellular CDMA pilot code on Q channel
    CDMA_UHF_SYNC_I        = 7,  // Cellular CDMA sync code on I channel
                                 //  (see N1 below)
    CDMA_UHF_SYNC_Q        = 8,  // Cellular CDMA sync code on Q channel
                                 //  (see N1 below)
    CDMA_UHF_PILOT_IQ      = 9,  // Cellular CDMA pilot code on I+Q channel
    CDMA_UHF_PILOT_SYNC_I  = 10, // Cellular CDMA pilot+sync combined code on
                                 //  I channel (see N1 below)
    CDMA_UHF_PILOT_SYNC_Q  = 11, // Cellular CDMA pilot+sync combined code on
                                 //  Q channel (see N1 below)
    CDMA_UHF_PILOT_SYNC_IQ = 12, // Cellular CDMA pilot+sync combined code on
                                 //  I+Q channel (see N1 below)
    SBAS_L1_I              = 13, // SBAS L1 on I channel
    // 14 unallocated to preserve legacy galileo values
    // The next new GenericType should use value 14.
    GALILEO_E1_BC          = 15, // Galileo E1 code (sum of E1B and E1C)
    GALILEO_E1_B           = 16, // E1B code for Galileo E1
    GALILEO_E1_C           = 17, // E1C code for Galileo E1
    NUM_GENERIC_TYPES      = 18,
    UNDEFINED_GENERIC_TYPE = 19
  };

  enum Group {
    PRIMARY                = 0,
    ALT1                   = 1,
    NUM_GROUPS             = 2,
    UNDEFINED_GROUP        = 3
    // add ALT2, ALT3, etc. as needed, up to a maximum of ALT7
  };

  enum DataModulationType {
    PILOT,
    DATA,
    PILOT_AND_DATA
  };

  enum Frequency {
    L1,
    L2,
    UHF,
    NUM_FREQUENCIES,
    UNDEFINED_FREQUENCY
  };

  SignalType()
    : SignalType(UNDEFINED_GENERIC_TYPE, UNDEFINED_GROUP) {}
  SignalType(GenericType genericType)
    : SignalType(genericType, PRIMARY) {}
  SignalType(GenericType genericType, Group group)
    : encodedValue_((static_cast<u8>(group) << GROUP_START_POS)
                   | static_cast<u8>(genericType)) {}
  SignalType(u8 serializedValue)
    : encodedValue_(serializedValue) {}
  SignalType(const std::string& str);
  static constexpr int NUM_SIGNAL_TYPES() {
    return static_cast<int>(
      (static_cast<u8>(NUM_GROUPS) << GROUP_START_POS)
      | static_cast<u8>(NUM_GENERIC_TYPES));
  }
  bool defined() const;
  GenericType genericType() const {
    return static_cast<GenericType>(encodedValue_ & GENERIC_TYPE_MASK);
  }
  Group group() const {
    return static_cast<Group>(
      (encodedValue_ & GROUP_MASK) >> GROUP_START_POS);
  }
  System system() const;
  bool isCombinationValid(const TxId& txId) const;
  bool isPrimitive() const;
  bool isTrackable() const;
  bool isAlt() const;
  bool isSbas() const;
  Frequency frequency() const;
  std::string toString() const;
  std::string toStringGroup() const;
  std::string toStringGenericType() const;
  u8 encodedValue() const { return encodedValue_;}
  // Overloading the () operator functions as a target-sensitive cast
  operator u8() const { return encodedValue_; }
  operator u32() const { return static_cast<u32>(encodedValue_); }
  operator s32() const { return static_cast<s32>(encodedValue_); }
  operator s64() const { return static_cast<s64>(encodedValue_); }
  inline bool operator==(const SignalType& rhs) const {
    return encodedValue_ == rhs.encodedValue_;
  }
  inline bool operator!=(const SignalType& rhs) const {
    return !operator==(rhs);
  }
  inline bool operator<(const SignalType& rhs) const {
    return encodedValue_ < rhs.encodedValue_;
  }

 private:
  static constexpr u8 GROUP_START_POS = 5;
  static constexpr u8 GROUP_MASK = 0b111 << GROUP_START_POS;
  static constexpr u8 GENERIC_TYPE_MASK = 0b11111;

  // Store Group (bits marked 'g') value and GenericType (bits marked 'T')
  // value inside encodedValue_:
  //   [ g g g T T T T T ]
  // Access to encodedValue_ should use GROUP_MASK, GROUP_START_POS, and
  // GENERIC_TYPE_MASK to facilitate the above bitpacking.
  //
  // encodedValue_ shall be the only data member of SignalType, so as to ensure
  // sizeof(SignalType)==1.
  u8 encodedValue_;
};

#endif

// Explanatory Notes:
//
// N1: As currently configured, accumulations from the CDMA SYNC channel are
// not used for tracking, but they are output to the output log files.
