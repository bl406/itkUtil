#pragma once

#include <itkRegionOfInterestImageFilter.h>
#include "itkFlatStructuringElement.h"
#include <itkConstantPadImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkChangeLabelImageFilter.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkNiftiImageIO.h>
#include <itkImportImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkResampleImageFilter.h>

// ReadImage
// WriteImage
// ImportImage
// BinaryThreshold
// ChangeLabel
// CropImage
// AutoCropImage
// ExtractImageROI
// ComputeDistanceMap
// ComputeLabelBoundingBox
// ConstantPad
// CastImageFilter
// ResampleImage
// FlipImage
namespace itkUtil
{ 
	template <typename InputImage, typename OutputImage>
	typename OutputImage::Pointer CastImage(const typename InputImage::Pointer &image)
	{
		typedef itk::CastImageFilter<InputImage, OutputImage> CastFitlerType;
		CastFitlerType::Pointer castFilter = CastFitlerType::New();
		castFilter->SetInput(image);
		castFilter->Update();
		return castFilter->GetOutput();
	}

	template <typename TImage>
	typename TImage::Pointer ConstantPad(const typename TImage::Pointer &image, const typename TImage::SizeType &lowerBound, const typename TImage::SizeType &upperBound,
		typename TImage::PixelType constant = itk::NumericTraits<typename TImage::PixelType>::Zero)
	{
		typedef itk::ConstantPadImageFilter<TImage, TImage> ConstantPadFilterType;
		ConstantPadFilterType::Pointer constantPad = ConstantPadFilterType::New();
		constantPad->SetInput(image);
		constantPad->SetPadLowerBound(lowerBound);
		constantPad->SetPadUpperBound(upperBound);
		constantPad->SetConstant(constant);
		constantPad->Update();

		return constantPad->GetOutput();
	}

	template <class TImage>
	typename TImage::Pointer ImportImage(const typename TImage::PixelType *pointer, const int *dim, const float *spacing, const float *origin)
	{
		const unsigned int dimension = TImage::ImageDimension;

		typename TImage::SpacingType imgSpacing;
		typename TImage::PointType imgOrigin;
		typename TImage::SizeType imgSize;
		typename TImage::RegionType imgRegion;
		for (unsigned int i = 0; i < dimension; ++i)
		{
			imgSpacing[i] = spacing[i];
			imgOrigin[i] = origin[i];
			imgSize[i] = dim[i];
		}
		imgRegion.SetSize(imgSize);

		typedef itk::ImportImageFilter<TImage::PixelType, TImage::ImageDimension> ImportFilterType;
		ImportFilterType::Pointer import = ImportFilterType::New();
		import->SetSpacing(imgSpacing);
		import->SetOrigin(imgOrigin);
		import->SetRegion(imgRegion);
		import->SetImportPointer((TImage::PixelType*)pointer, imgRegion.GetNumberOfPixels(), false);
		import->Update();

		typename TImage::Pointer image = import->GetOutput();
		image->DisconnectPipeline();
		return image;
	}
	
	template <typename TImage>
	typename TImage::Pointer ReadImage(const std::string fileName, const bool zeroOrigin = false)
	{
		typedef itk::ImageFileReader<TImage> ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(fileName.c_str());
		try
		{
			reader->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cout << "Caught an exception: " << std::endl;
			std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
			throw err;
		}
		catch (...)
		{
			std::cout << "Error while reading in image for patient " << fileName << std::endl;
			throw;
		}

		typename TImage::Pointer image = reader->GetOutput();
		if (zeroOrigin)
		{
			double origin[TImage::ImageDimension];
			for (unsigned int i = 0; i < TImage::ImageDimension; i++)
			{
				origin[i] = 0;
			}
			image->SetOrigin(origin);
		}
		return image;
	}
	
	template <class ImageType>
	void WriteImage(const typename ImageType::Pointer &image, const std::string &filename, bool compression = false)
	{

		typedef itk::ImageFileWriter< ImageType > WriterType;
		typename  WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(filename.c_str());
		writer->SetInput(image);
		writer->SetUseCompression(compression);

		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject &err)
		{
			std::cout << "Exception Object caught: " << std::endl;
			std::cout << err << std::endl;
			throw(err);
		}
	}

#if _MSC_VER >= 1600 
	template <typename TInputImage, typename TOutputImage = TInputImage>
	typename TOutputImage::Pointer ChangeLabel(const typename TInputImage::Pointer inputImage, 
									 const typename TInputImage::PixelType &original, const typename TOutputImage::PixelType &result)
#else
	template <typename TInputImage, typename TOutputImage>
	typename TOutputImage::Pointer ChangeLabel(const typename TInputImage::Pointer inputImage, 
									 const typename TInputImage::PixelType &original, const typename TOutputImage::PixelType &result)
#endif
	{
		typedef itk::ChangeLabelImageFilter<TInputImage, TOutputImage> ChangeLabelFilterType;
		ChangeLabelFilterType::Pointer changeLabel = ChangeLabelFilterType::New();
		changeLabel->SetChange(original, result);
		changeLabel->SetInput(inputImage);
		changeLabel->Update();
		return changeLabel->GetOutput();
	}

#if _MSC_VER >= 1600  
	template <typename TInputImage, typename TOutputImage = TInputImage>
	typename TOutputImage::Pointer CropImage(const typename TInputImage::Pointer inputImage,
		const typename TInputImage::RegionType& roiRegion)
#else
	template <typename TInputImage, typename TOutputImage>
	typename TOutputImage::Pointer CropImage(const typename TInputImage::Pointer inputImage,
		const typename TInputImage::RegionType& roiRegion)
#endif
	{
		typedef itk::ExtractImageFilter<TInputImage, TOutputImage> ExtractFilterType;
		ExtractFilterType::Pointer extract = ExtractFilterType::New();
		extract->SetInput(inputImage);
		extract->SetExtractionRegion(roiRegion);
		extract->Update();

		TOutputImage::Pointer roi = extract->GetOutput();
		roi->DisconnectPipeline();
		return roi;
	}

	template <class TImage,class TLabelImage>
	typename TImage::Pointer AutoCropImage(const typename TImage::Pointer inputImage,
										  const typename TLabelImage::Pointer labelImage,
										  typename TLabelImage::PixelType label,
										  const typename TLabelImage::SizeType &padding,
										  bool shiftorigin = true)
	{
		return CropImage<ImageType>(inputImage, ComputeLabelBound<TLabelImage>(labelImage, label, padding), shiftorigin);
	}
	
	template <class TImage>
	typename TImage::Pointer ExtractImageROI(const typename TImage::Pointer input,
		const typename TImage::RegionType &roi)
	{
		typedef itk::RegionOfInterestImageFilter<TImage, TImage> ROIFilterType;
		ROIFilterType::Pointer roiFilter = ROIFilterType::New();
		roiFilter->SetInput(input);
		roiFilter->SetRegionOfInterest(roi);

		try
		{
			roiFilter->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cout << "Caught an exception: " << std::endl;
			std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
			throw err;
		}
		catch (...)
		{
			std::cout << "ExtractImageROI encount an error! " << std::endl;
			throw;
		}

		TImage::Pointer output = roiFilter->GetOutput();
		output->DisconnectPipeline();
		return output;
	}


	template<typename LabelVolumeType>
	itk::Image<float, 3>::Pointer ComputeDistanceMap(typename LabelVolumeType::Pointer lv, bool insidePositive = true)
	{
		typedef itk::Image<float, 3> FloatVolumeType;
		itk::SignedMaurerDistanceMapImageFilter<LabelVolumeType, FloatVolumeType>::Pointer distanceFilter = itk::SignedMaurerDistanceMapImageFilter<LabelVolumeType, FloatVolumeType>::New();
		distanceFilter->SetInput(lv);
		distanceFilter->SetInsideIsPositive(insidePositive);
		distanceFilter->Update();
		return distanceFilter->GetOutput();
	}

	template <class TLabelImage>
	typename TLabelImage::RegionType ComputeLabelBoundingBox(const typename TLabelImage::Pointer labelImage,
		typename TLabelImage::PixelType label,
		const typename TLabelImage::SizeType &padding)
	{
		typedef itk::ImageRegionIteratorWithIndex<TLabelImage> ItType;
		ItType it(labelImage, labelImage->GetBufferedRegion());
		TLabelImage::IndexType minid, maxid;
		maxid = it.GetRegion().GetIndex();
		minid = it.GetRegion().GetUpperIndex();
	
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (it.Get() == label)
			{
				for (int i = 0; i < 3; ++i)
				{
					if (it.GetIndex()[i] > maxid[i]) maxid[i] = it.GetIndex()[i];
					if (it.GetIndex()[i] < minid[i]) minid[i] = it.GetIndex()[i];
				}
			}
		}

		TLabelImage::RegionType roiRegion;
		minid -= padding;
		maxid += padding;
		for (int i = 0; i < 3; ++i)
		{
			if (minid[i]<it.GetRegion().GetIndex()[i]) minid[i] = it.GetRegion().GetIndex()[i];
			if (maxid[i] >= it.GetRegion().GetIndex()[i] + it.GetRegion().GetSize()[i])
				maxid[i] = it.GetRegion().GetIndex()[i] + it.GetRegion().GetSize()[i] - 1;
		}
		roiRegion.SetIndex(minid);
		roiRegion.SetUpperIndex(maxid);

		return roiRegion;
	}

#if _MSC_VER >= 1600 
	template<typename PixelType, int InputDimension, int OutputDiemnsion = InputDimension>
	static typename itk::Image<PixelType, OutputDiemnsion>::Pointer
		ExtractImageFromVolume(typename itk::Image<PixelType, InputDimension>::Pointer inputVolume, const itk::ImageRegion<InputDimension> &region)
#else 
	template<typename PixelType, int InputDimension, int OutputDiemnsion>
	static typename itk::Image<PixelType, OutputDiemnsion>::Pointer
		ExtractImageFromVolume(typename itk::Image<PixelType, InputDimension>::Pointer inputVolume, const itk::ImageRegion<InputDimension> &region)

#endif
	{
		Q_ASSERT(OutputDiemnsion <= InputDimension);
		typedef itk::Image<PixelType, InputDimension> InputImageType;
		typedef itk::Image<PixelType, OutputDiemnsion> OutputImageType;
		typedef itk::ExtractImageFilter<InputImageType, OutputImageType> ExtractImageFilterType;
		ExtractImageFilterType::Pointer extractor = ExtractImageFilterType::New();
		extractor->SetInput(inputVolume);
		extractor->SetExtractionRegion(region);
		extractor->SetDirectionCollapseToIdentity();
		extractor->Update();

		return extractor->GetOutput();
	}


	template <class TImage, class TTransform, class TReferImage>
	typename TImage::Pointer ResampleImage(const typename TImage::Pointer input,
		const typename TReferImage::Pointer refer, 
		const typename TTransform::Pointer transform)
	{
		typedef itk::ResampleImageFilter<TImage,TImage> ResampleFilterType;
		ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();

		
		resampleFilter->SetTransform(transform);
		resampleFilter->SetInput(input);
		resampleFilter->SetSize(refer->GetLargestPossibleRegion().GetSize());		
		resampleFilter->SetOutputStartIndex(refer->GetLargestPossibleRegion().GetIndex());
		resampleFilter->SetOutputOrigin(refer->GetOrigin());
		resampleFilter->SetOutputSpacing(refer->GetSpacing());
		resampleFilter->SetOutputDirection(refer->GetDirection());

		try
		{
			resampleFilter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cout << "Caught an exception: " << std::endl;
			std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
			throw err;
		}
		catch(...)
		{
			std::cout << "ResampleImage encount an error! " << std::endl;
			throw;
		}
	}
	
		// Not like the itkFlipImageFilter, FlipImage flips the data memory.
	template <class TImage>
	typename TImage::Pointer FlipImage(const typename TImage::Pointer input, 
				const typename TImage::IndexType &flipAxes, 
				bool flipmemory = true, bool flipaboutorigin = false)
	{
		typedef itk::FlipImageFilter<TImage> FlipImageFilter;
		FlipImageFilter::Pointer flip = FlipImageFilter::New();
		flip->SetInput(input);
		
		FlipImageFilter::FlipAxesArrayType axesArray;
		axesArray.Fill(0);
		for(unsigned int i = 0; i < TImage::ImageDimension; ++i){
			axesArray[i] = flipAxes[i];			
		}

		flip->SetFlipAxes(axesArray);
		flip->SetFlipAboutOrigin(flipaboutorigin);
		flip->Update();

		if(flipmemory)
		{
			typedef itk::ResampleImageFilter<TImage, TImage, float> ResampleFilterType;
			ResampleFilterType::Pointer resample = ResampleFilterType::New();
			resample->SetInput(input);
			resample->SetReferenceImage(flip->GetOutput());
			resample->SetUseReferenceImage(true);
			resample->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<TImage,float>::New());
			resample->Update();

			TImage::Pointer output = resample->GetOutput();
			output->DisconnectPipeline();

			output->SetRegions(input->GetBufferedRegion());
			output->SetOrigin(input->GetOrigin());
			output->SetDirection(input->GetDirection());
			
			return output;
		}
		else
		{
			return flip->GetOutput();
		}
	}

	class ShowProgressObject
	{
	public:
		ShowProgressObject(itk::ProcessObject* o)
		{
			m_Process = o;
		}
		void ShowProgress()
		{
			std::cout << "Progress " << m_Process->GetProgress() << std::endl;
		}
		itk::ProcessObject::Pointer m_Process;
	};
};