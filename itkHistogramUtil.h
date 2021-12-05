#pragma once

#include <itkScalarImageToHistogramGenerator.h> 
#include <itkImage.h>
#include <itkHistogram.h>

namespace itkUtil
{
	typedef itk::Statistics::Histogram<double> HistogramType;
	typedef HistogramType::Pointer HistogramPointer;

	template <typename InputImage>
	const HistogramType* ComputeHistogram(typename InputImage::Pointer input, int nbin = 256)
	{
		typedef itk::Statistics::ScalarImageToHistogramGenerator<InputImage> HistogramGenerType;
		HistogramGenerType::Pointer histGener = HistogramGenerType::New();
		histGener->SetInput(input);
		histGener->SetNumberOfBins(nbin);
		histGener->Compute();

		return histGener->GetOutput();
	}

	template <typename InputImage>
	typename InputImage::PixelType ComputeIntensityQuantile(typename InputImage::Pointer input, double quantile, int nbin = 256)
	{
		typedef itk::Statistics::ScalarImageToHistogramGenerator<InputImage> HistogramGenerType;
		HistogramGenerType::Pointer histGener = HistogramGenerType::New();
		histGener->SetInput(input);
		histGener->SetNumberOfBins(nbin);
		histGener->Compute();

		return histGener->GetOutput()->Quantile(0, quantile);
	}
}