#pragma once

#include "itkUtil.h"
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBlackTopHatImageFilter.h>

namespace itkUtil{
	inline void fixRadius(int &xrad, int &yrad, int &zrad)
	{
		if (yrad < 0)
			yrad = xrad;
		if (zrad < 0)
			zrad = xrad;
	}
	
	template <class TImage>
	typename TImage::Pointer BinaryOpening(const typename TImage::Pointer input, typename TImage::PixelType foreground, int xrad,
		int yrad = -1, int zrad = -1)
	{
		const unsigned int dim = TImage::ImageDimension;
		typedef typename itk::FlatStructuringElement< dim > SRType;

		fixRadius(xrad, yrad, zrad);

		typename SRType::RadiusType rad;
		rad[0] = xrad;
		rad[1] = yrad;
		if (dim > 2)
			rad[2] = zrad;

		SRType kernel;

		kernel = SRType::Ball(rad);

		typedef typename itk::BinaryMorphologicalOpeningImageFilter<TImage, TImage, SRType> FiltType;
		typename FiltType::Pointer filter = FiltType::New();
		filter->SetInput(input);
		filter->SetKernel(kernel);
		filter->SetForegroundValue(foreground);

		typename TImage::Pointer result = filter->GetOutput();
		result->Update();
		result->DisconnectPipeline();
		return(result);
	}
	
	template <class TImage>
	typename TImage::Pointer BinaryClosing(const typename TImage::Pointer input, typename TImage::PixelType foreground, int xrad,
		int yrad = -1, int zrad = -1)
	{
		const unsigned int dim = TImage::ImageDimension;
		typedef typename itk::FlatStructuringElement< dim > SRType;

		fixRadius(xrad, yrad, zrad);

		typename SRType::RadiusType rad;
		rad[0] = xrad;
		rad[1] = yrad;
		if (dim > 2)
			rad[2] = zrad;

		SRType kernel;

		kernel = SRType::Ball(rad);

		typedef typename itk::BinaryMorphologicalClosingImageFilter<TImage, TImage, SRType> FiltType;
		typename FiltType::Pointer filter = FiltType::New();
		filter->SetInput(input);
		filter->SetKernel(kernel);
		filter->SetForegroundValue(foreground);

		typename TImage::Pointer result = filter->GetOutput();
		result->Update();
		result->DisconnectPipeline();
		return(result);
	}

	template <class TImage>
	typename TImage::Pointer BlackTopHat(const typename TImage::Pointer input, int xrad, int yrad = -1, int zrad = -1)
	{
		const unsigned int dim = TImage::ImageDimension;
		typedef typename itk::FlatStructuringElement< dim > SRType;

		fixRadius(xrad, yrad, zrad);

		typename SRType::RadiusType rad;
		rad[0] = xrad;
		rad[1] = yrad;
		if (dim > 2)
			rad[2] = zrad;

		SRType kernel;

		kernel = SRType::Ball(rad);

		typedef typename itk::BlackTopHatImageFilter<TImage, TImage, SRType> FiltType;
		typename FiltType::Pointer filter = FiltType::New();
		filter->SetInput(input);
		filter->SetKernel(kernel);

		typename TImage::Pointer result = filter->GetOutput();
		result->Update();
		result->DisconnectPipeline();
		return(result);
	}
};