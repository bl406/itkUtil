#pragma once

#include <itkBinaryThresholdImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>

namespace itkUtil
{
	template <typename TInputImage, typename TOutputImage>
	typename TOutputImage::Pointer BinaryThreshold(typename TInputImage::Pointer input,
		typename TInputImage::PixelType lowThresh, typename TInputImage::PixelType upperThresh,
		typename TOutputImage::PixelType insideValue = static_cast<TOutputImage::PixelType>(1),
		typename TOutputImage::PixelType outsideValue = static_cast<TOutputImage::PixelType>(0))
	{
		typedef itk::BinaryThresholdImageFilter<TInputImage, TOutputImage> ThreshodFilterType;
		ThreshodFilterType::Pointer threshold = ThreshodFilterType::New();
		threshold->SetInput(input);
		threshold->SetLowerThreshold(lowThresh);
		threshold->SetUpperThreshold(upperThresh);
		threshold->SetInsideValue(insideValue);
		threshold->SetOutsideValue(outsideValue);
		threshold->Update();
		return threshold->GetOutput();
	}

	template <typename TInputImage, typename TOutputImage>
	typename TOutputImage::Pointer ConnectedThreshold(typename TInputImage::Pointer input, const typename TInputImage::IndexType &seed,
		typename TInputImage::PixelType lowThresh, typename TInputImage::PixelType upperThresh,
		typename TOutputImage::PixelType replaceValue = static_cast<TOutputImage::PixelType>(1))
	{
		typedef itk::ConnectedThresholdImageFilter<TInputImage, TOutputImage>  ConnectedThresholdFilterType;
		ConnectedThresholdFilterType::Pointer filter = ConnectedThresholdFilterType::New();
		filter->SetInput(input);
		filter->SetSeed(seed);
		filter->SetLower(lowThresh);
		filter->SetUpper(upperThresh);
		filter->SetReplaceValue(replaceValue);
		filter->Update();
		return filter->GetOutput();
	}
};