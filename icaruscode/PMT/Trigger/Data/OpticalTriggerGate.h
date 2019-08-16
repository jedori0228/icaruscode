/**
 * @file   icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h
 * @brief  A trigger gate data object for optical detector electronics.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_OPTICALTRIGGERGATE_H
#define ICARUSCODE_PMT_TRIGGER_DATA_OPTICALTRIGGERGATE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // detinfo::timescales
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <utility> // std::move()


//------------------------------------------------------------------------------

//
// declarations
//
namespace icarus::trigger {
  class OpticalTriggerGate;
  std::ostream& operator<< (std::ostream&, OpticalTriggerGate const&);
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Logical multi-level gate associated to one or more waveforms.
 * 
 * This object is a trigger gate associated with one or more optical waveforms.
 */
class icarus::trigger::OpticalTriggerGate
  : public icarus::trigger::TriggerGateData
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_time_ticks>
{
  
    public:
  /// Type for gate data access.
  using GateData_t = icarus::trigger::TriggerGateData
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_time_ticks>
    ;
  
  
  /// Constructor: a closed gate with no associated waveform (`add()` them).
  OpticalTriggerGate() = default;
  
  /// Constructor: a closed gate for the channel in `waveform`.
  OpticalTriggerGate(raw::OpDetWaveform const& waveform)
    : GateData_t()
    , fWaveforms({ &waveform })
    {}
  
  /// Adds another waveform to the gate (unless it has already been added).
  bool add(raw::OpDetWaveform const& waveform);
  
  // --- BEGIN Query -----------------------------------------------------------
  /// @name Query
  /// @{
  
  /// Access to the underlying gate level data (mutable).
  GateData_t& gateLevels() { return *this; }
  
  /// Access to the underlying gate level data (immutable).
  GateData_t const& gateLevels() const { return *this; }
  
  /// Returns a list of pointers to the waveforms associated to the gate,
  /// sorted.
  std::vector<raw::OpDetWaveform const*> waveforms() const
    { return fWaveforms; }
  
  // --- END Query -------------------------------------------------------------
  
  
  // --- BEGIN Combination operations ------------------------------------------
  /// @name Combination operations
  /// @{
  
  /**
   * @brief Combines with a gate, keeping the minimum opening among the two.
   * @param other gate to combine to
   * @return this object
   * @see `Max()`
   * 
   * Multi-level equivalent of an _and_ logical operation.
   */
  OpticalTriggerGate& Min(OpticalTriggerGate const& other);

  /**
   * @brief Combines with a gate, keeping the maximum opening among the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Sum()`
   * 
   * Multi-level equivalent of an _or_ logical operation.
   */
  OpticalTriggerGate& Max(OpticalTriggerGate const& other);

  /**
   * @brief Combines with a gate, keeping the sum of openings of the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Max()`
   */
  OpticalTriggerGate& Sum(OpticalTriggerGate const& other);

  /**
   * @brief Combines with a gate, keeping the product of openings of the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Max()`, `Sum()`
   */
  OpticalTriggerGate& Mul(OpticalTriggerGate const& other);

  /**
   * @brief Returns a gate with the minimum opening between the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the minimum opening among `a` and `b`
   * @see `Max()`
   * 
   * Multi-level equivalent of an _and_ logical operation.
   */
  static OpticalTriggerGate Min
    (OpticalTriggerGate const& a, OpticalTriggerGate const& b);

  /**
   * @brief Returns a gate with the maximum opening between the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the maximum opening among `a` and `b`
   * @see `Min()`, `Sum()`
   * 
   * Multi-level equivalent of an _or_ logical operation.
   */
  static OpticalTriggerGate Max
    (OpticalTriggerGate const& a, OpticalTriggerGate const& b);

  /**
   * @brief Returns a gate with opening sum of the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the total opening of `a` and `b`
   * @see `Max()`
   */
  static OpticalTriggerGate Sum
    (OpticalTriggerGate const& a, OpticalTriggerGate const& b);

  /**
   * @brief Returns a gate with opening product of the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the product of openings of `a` and `b`
   * @see `Max()`
   */
  static OpticalTriggerGate Mul
    (OpticalTriggerGate const& a, OpticalTriggerGate const& b);


  /**
   * @brief Returns a gate combination of the openings of two other gates.
   * @tparam Op binary operation: `OpeningCount_t` (x2) to `OpeningCount_t`
   * @param op symmetric binary combination operation
   * @param a first gate
   * @param b second gate
   * @return gate with opening combination of `a` and `b`
   * 
   * For this algorithm to work, the operation needs to be symmetric, i.e.
   * `op(c1, c2) == op(c2, c1)` for every valid combinations of counts
   * `c1` and `c2`.
   * 
   */
  template <typename Op>
  static OpticalTriggerGate SymmetricCombination(
    Op&& op, OpticalTriggerGate const& a, OpticalTriggerGate const& b,
    ClockTicks_t aDelay = ClockTicks_t{ 0 },
    ClockTicks_t bDelay = ClockTicks_t{ 0 }
    );
  
  /// @}
  // --- END Combination operations --------------------------------------------
  
  
    protected:
  // we allow some manipulation by the derived classes
  
  /// Internal list of registered waveforms.
  using Waveforms_t = std::vector<raw::OpDetWaveform const*>;
  
  /// Protected constructor: set the data directly.
  OpticalTriggerGate(GateEvolution_t&& gateLevel, Waveforms_t&& waveforms)
    : TriggerGateData(std::move(gateLevel)), fWaveforms(std::move(waveforms))
    {}
  
  
  /// Registers the waveforms from the specified list.
  void registerWaveforms(Waveforms_t const& moreWaveforms);
  
  /// Registers the waveforms from the `other` gate into this one.
  void mergeWaveformsFromGate(OpticalTriggerGate const& other)
    { registerWaveforms(other.waveforms()); }
  
  
  /// Registers the waveforms from the `other` gate into this one.
  static Waveforms_t mergeWaveforms(Waveforms_t const& a, Waveforms_t const& b);
  
    private:
  
  /// List of waveforms involved in this channel.
  std::vector<raw::OpDetWaveform const*> fWaveforms;
  
  
}; // class icarus::trigger::OpticalTriggerGate


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_DATA_OPTICALTRIGGERGATE_H