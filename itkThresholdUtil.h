#pragma once

#include "itkUtil.h"
#include <itkConnectedThresholdImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>

namespace itkUtil{
	
	template <typename TInputImage, typename TOutputImage>
	typename TOutputImage::Pointer OtsuThreshold(
		typename TInputImage::Pointer input,
		typename TOutputImage::PixelType insideValue = static_cast<TOutputImage::PixelType>(1),
		typename TOutputImage::PixelType outsideValue = static_cast<TOutputImage::PixelType>(0)
		)
	{
		typedef itk::OtsuThresholdImageFilter<TInputImage, TOutputImage> OtsuThresholdFilterType;
		OtsuThresholdFilterType::Pointer otsu = OtsuThresholdFilterType::New();
		otsu->SetInput(input);
		otsu->SetInsideValue(insideValue);
		otsu->SetOutsideValue(outsideValue);
		otsu->Update();
		
		return otsu->GetOutput();
	}
	
	template <typename TInputImage, typename TOutputImage>
	typename TOutputImage::Pointer ConnectedThreshold(
		typename TInputImage::Pointer input,
		typename TInputImage::PixelType lowThresh, 
		typename TInputImage::PixelType upperThresh,
		typename TOutputImage::PixelType insideValue = static_cast<TOutputImage::PixelType>(1),
		typename TOutputImage::PixelType outsideValue = static_cast<TOutputImage::PixelType>(0)
		)
	{
		typedef itk::BinaryThresholdImageFilter<TInputImage, TOutputImage> ThreshodFilterType;
		ThreshodFilterType::Pointer threshold = ThreshodFilterType::New();
		threshold->SetInput(input);
		threshold->SetLowerThreshold(upperThresh);
		threshold->SetUpperThreshold(std::numeric_limits<TInputImage::PixelType>::max());
		threshold->SetInsideValue(insideValue);
		threshold->SetOutsideValue(outsideValue);
		threshold->Update();

		typedef itk::ConnectedThresholdImageFilter<TInputImage, TOutputImage> ConnectedThresholdType;
		ConnectedThresholdType::Pointer connectedThreshold = ConnectedThresholdType::New();
		connectedThreshold->SetInput(input);
		connectedThreshold->SetUpper(std::numeric_limits<TInputImage::PixelType>::max());
		connectedThreshold->SetLower(lowThresh);
		itk::ImageRegionConstIteratorWithIndex<TOutputImage> it(threshold->GetOutput(), threshold->GetOutput()->GetBufferedRegion());
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (it.Get() > 0) connectedThreshold->AddSeed(it.GetIndex());
		}
		connectedThreshold->Update();
		return connectedThreshold->GetOutput();
	}
	
};