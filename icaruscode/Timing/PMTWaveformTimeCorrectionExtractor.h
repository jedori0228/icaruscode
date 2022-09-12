/**
 * @file   TimingCorrectionExtraction_module.cc
 * @brief  Extract timing correction and adjust waveform starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

#ifndef ICARUSCODE_TIMING_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H
#define ICARUSCODE_TIMING_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H

// ICARUS/SBN libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/Timing/DataProducts/PMTWaveformTimeCorrection.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

namespace icarus::timing { class PMTWaveformTimeCorrectionExtractor; }

class icarus::timing::PMTWaveformTimeCorrectionExtractor {
	
	public: 

        PMTWaveformTimeCorrectionExtractor(
            detinfo::DetectorClocksData const detTimingService,
            icarusDB::IICARUSChannelMap const & channelMapService,
            icarusDB::PMTTimingCorrections & pmtTimingCorrectionsService, 
            bool const & verbose );

        ~PMTWaveformTimeCorrectionExtractor(){};

        void findWaveformTimeCorrections(   
            raw::OpDetWaveform const & wave,
            std::string const & cateogry, 
            unsigned int const & waveChannelID, 
            bool const & correctCableDelay,
            std::vector<PMTWaveformTimeCorrection> & corrections );


	private:

        detinfo::DetectorClocksData const fClocksData;

        icarusDB::IICARUSChannelMap const & fChannelMap;

        icarusDB::PMTTimingCorrections & fPMTTimingCorrectionsService;

        bool const & fVerbose;

        std::map<unsigned int, std::vector<unsigned int>> 
            fCrateFragmentMap {
                {0x0070 , { 0 , 1 , 2  }}, 
                {0x0060 , { 3 , 4 , 5  }}, 
                {0x0050 , { 6 , 7 , 8  }}, 
                {0x0040 , { 9 , 10, 11 }},
                {0x0030 , { 12, 13, 14 }}, 
                {0x0020 , { 15, 16, 17 }}, 
                {0x0010 , { 18, 19, 20 }}, 
                {0x0000 , { 21, 22, 23 }}, 
            }; 

        template<typename T>
            size_t getMaxBin( 
                std::vector<T> vv, 
                size_t startElement, 
                size_t endElement);

        template<typename T>
            size_t getMinBin( 
                std::vector<T> vv, 
                size_t startElement, 
                size_t endElement);

        template<typename T>
            size_t getStartSample( std::vector<T> vv );

};

#endif //ICARUSCODE_TIMING_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H