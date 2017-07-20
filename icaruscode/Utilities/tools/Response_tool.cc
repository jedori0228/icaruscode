////////////////////////////////////////////////////////////////////////
/// \file   Response.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/SignalShaping.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Utilities/make_tool.h"
#include "icaruscode/Utilities/tools/IWaveformTool.h"
#include "IFieldResponse.h"
#include "IElectronicsResponse.h"
#include "IFilter.h"

#include <fstream>
#include <iomanip>

namespace icarus_tool
{

class Response : IResponse
{
public:
    explicit Response(const fhicl::ParameterSet& pset);
    
    ~Response() {}
    
    void configure(const fhicl::ParameterSet& pset)   override;
    void setResponse(double weight)                   override;
    void outputHistograms(art::TFileDirectory&) const override;
    
    size_t                      getPlane()               const override {return fThisPlane;}
    
    const IFieldResponse*       getFieldResponse()       const override {return fFieldResponse.get();}
    const IElectronicsResponse* getElectronicsResponse() const override {return fElectronicsResponse.get();}
    const IFilter*              getFilter()              const override {return fFilter.get();}
    
    const util::SignalShaping&  getSignalShaping()       const override {return fSignalShaping;}
    
private:
    using IFieldResponsePtr       = std::unique_ptr<icarus_tool::IFieldResponse>;
    using IElectronicsResponsePtr = std::unique_ptr<icarus_tool::IElectronicsResponse>;
    using IFilterPtr              = std::unique_ptr<icarus_tool::IFilter>;

    // Utility routine for converting numbers to strings
    std::string             numberToString(int number);
    
    // Member variables from the fhicl file
    size_t                  fThisPlane;
    double                  f3DCorrection;
    double                  fTimeScaleFactor;
    int                     fDeconvPol;
    
    // Keep track of our base tools
    IFieldResponsePtr       fFieldResponse;
    IElectronicsResponsePtr fElectronicsResponse;
    IFilterPtr              fFilter;
    
    // The actual response function
    util::SignalShaping     fSignalShaping;
};
    
//----------------------------------------------------------------------
// Constructor.
Response::Response(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
void Response::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fThisPlane       = pset.get<size_t>("Plane");
    f3DCorrection    = pset.get<size_t>("Correction3D");
    fTimeScaleFactor = pset.get<size_t>("TimeScaleFactor");
    fDeconvPol       = pset.get<int>("DeconvPol");
    
    // Build out the underlying tools we'll be using
    fFieldResponse       = art::make_tool<icarus_tool::IFieldResponse>(pset.get<fhicl::ParameterSet>("FieldResponse"));
    fElectronicsResponse = art::make_tool<icarus_tool::IElectronicsResponse>(pset.get<fhicl::ParameterSet>("ElectronicsResponse"));
    fFilter              = art::make_tool<icarus_tool::IFilter>(pset.get<fhicl::ParameterSet>("Filter"));
    
    return;
}
    
void Response::setResponse(double weight)
{
    // We'll need the FFT service
    art::ServiceHandle<util::LArFFT> fastFourierTransform;
    
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Recover the current set up
    std::string fftOptions  = fastFourierTransform->FFTOptions();
    size_t      nFFTFitBins = fastFourierTransform->FFTFitBins();
    size_t      fftSizeIn   = fastFourierTransform->FFTSize();
    size_t      fftSize     = fftSizeIn;
    
    // First of all set the field response
    fFieldResponse->setResponse(weight, f3DCorrection, fTimeScaleFactor);
    
    // Make sure the FFT can handle this
    size_t nFieldBins = fFieldResponse->getNumBins();
    
    // Reset the FFT if it is not big enough to handle current size
    if (nFieldBins * 4 > fftSize)
    {
        fftSize = 4 * nFieldBins;
        
        fastFourierTransform->ReinitializeFFT( fftSize, fftOptions, nFFTFitBins);
    }
        
    // handle the electronics response for this plane
    fElectronicsResponse->setResponse(fftSize, fFieldResponse->getBinWidth());
    
    // Set up the filter
    fFilter->setResponse(fftSizeIn, f3DCorrection, fTimeScaleFactor);
    
    // Add these elements to the SignalShaping class
    fSignalShaping.Reset();
    fSignalShaping.AddResponseFunction(fFieldResponse->getResponseVec());
    fSignalShaping.AddResponseFunction(fElectronicsResponse->getResponseVec());
    fSignalShaping.save_response();
    fSignalShaping.set_normflag(false);
    
    // Now set to the task of determing the actual sampling response
    // We hve to remember that the bin size for determining the field response probably
    // does not match that for the detector readout so we'll need to "convert"
    // from one to the other.
    std::vector<double> samplingTimeVec( fftSize, 0. );
    
    // Recover the combined response from above
    const std::vector<double>& curResponseVec = fSignalShaping.Response_save();
    
    // Need two factors: 1) the detector sampling rate and 2) the response sampling rate
    double samplingRate = detprop->SamplingRate() * 1.e-3;       // We want this in us/bin
    double responseRate = fFieldResponse->getBinWidth() * 1.e-3; // We want this in us/bin
    double rateRatio    = samplingRate / responseRate;
    
    // The idea is to step through each bin of the sampling response vector and then to
    // look up the corresponding bins in the current response vector. Since the two sample
    // rates are not the same there will be some "stretching" between the two. In addition,
    // we want to continue to allow for the possibility for further sample stretching
    double binScaleFactor = rateRatio * f3DCorrection * fTimeScaleFactor;
    
    // ok, do the loop
    for(size_t sampleIdx = 0; sampleIdx < samplingTimeVec.size(); sampleIdx++)
    {
        // calculate the index for the response
        size_t responseLowIdx = std::floor(sampleIdx * binScaleFactor);
        
        if (responseLowIdx < curResponseVec.size())
        {
            // Calculate the index for the next bin
            size_t responseHiIdx = std::floor((sampleIdx + 1) * binScaleFactor);
            
            // This can't happen? But protect against zero divides...
            if (responseHiIdx == responseLowIdx) responseHiIdx += 1;
            
            if (responseHiIdx < curResponseVec.size())
            {
                
                // Now interpolate between the two bins to get the sampling response for this bin
                double responseSlope = (curResponseVec.at(responseHiIdx) - curResponseVec.at(responseLowIdx)) / (responseHiIdx - responseLowIdx);
                double response      = curResponseVec.at(responseLowIdx) + 0.5 * responseSlope * (responseHiIdx - responseLowIdx);
                
                samplingTimeVec.at(sampleIdx) = response;
            }
        }
    }

    fSignalShaping.AddResponseFunction( samplingTimeVec, true);
    
    // Currently we only have fine binning "fFieldBinWidth"
    // for the field and electronic responses.
    // Now we are sampling the convoluted field-electronic response
    // with the nominal sampling.
    // We may consider to do the same for the filters as well.
    if (fftSizeIn != fftSize) fastFourierTransform->ReinitializeFFT(fftSizeIn, fftOptions, nFFTFitBins);
    
    // Finalize the Signal Shaping
    fSignalShaping.AddFilterFunction(fFilter->getResponseVec());
    fSignalShaping.SetDeconvKernelPolarity( fDeconvPol );
    fSignalShaping.CalculateDeconvKernel();
    
    return;
}
    
void Response::outputHistograms(art::TFileDirectory& histDir) const
{
    // Create a subfolder in which to place the "response" histograms
    std::string thisResponse = "ResponsesPlane_" + std::to_string(fThisPlane);
    
    art::TFileDirectory dir = histDir.mkdir(thisResponse.c_str());
    
    // Do the field response histograms
    fFieldResponse->outputHistograms(dir);
    fElectronicsResponse->outputHistograms(dir);
    fFilter->outputHistograms(dir);
    
    // Now make hists for the full response
    std::string dirName = "Response_" + std::to_string(fThisPlane);
    
    art::TFileDirectory        responesDir  = dir.mkdir(dirName.c_str());
    const std::vector<double>& responseVec  = this->getSignalShaping().Response();
    double                     numBins      = responseVec.size();
    std::string                histName     = "Response_Plane_" + std::to_string(fThisPlane);
    TH1D*                      hist         = dir.make<TH1D>(histName.c_str(), "Response;Time(ticks)", numBins, 0., numBins);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        hist->Fill(bin, responseVec.at(bin));
    }

    // Get the FFT, need the waveform tool for consistency
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    std::unique_ptr<icarus_tool::IWaveformTool> waveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);
    
    std::vector<double> powerVec;
    
    waveformTool->getFFTPower(responseVec, powerVec);
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      maxFreq      = 500. / samplingRate;
    double      freqWidth    = maxFreq / powerVec.size();
    std::string freqName     = "Response_FFTPlane_" + std::to_string(fThisPlane);
    TH1D*       freqHist     = dir.make<TH1D>(freqName.c_str(), "Response;Frequency(MHz)", powerVec.size(), 0., maxFreq);
    
    for(size_t idx = 0; idx < powerVec.size(); idx++)
    {
        double freq = freqWidth * (idx + 0.5);
        
        freqHist->Fill(freq, powerVec.at(idx));
    }
    

    return;
}

std::string Response::numberToString(int number)
{
    std::ostringstream string;
    
    string << std::setfill('0') << std::setw(2) << number;
    
    return string.str();
}

    
DEFINE_ART_CLASS_TOOL(Response)
}