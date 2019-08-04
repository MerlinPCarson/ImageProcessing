///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//						Modified:   Merlin Carson
//						Date:	    Fall 2018
//						Why:	    Programmed algorithms-
//								Uniform Quantization,
//								Populosity
//								Naive Threshold Dithering
//								Brightness Preserving Threshold
//								Random Dithering
//								Clustered Dithering
//								Floyd-Steinberg Dithering
//								Box Filter
//								Bartlett Filter
//								Gaussian Filter
//								Compositing Over
//								Compositing Inside
//								Compositing Outside
//								Compositing Atop
//								Compositing Xor
//															
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <algorithm> 

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
	cout << "Grayscale" << endl;
	double grayVal = 0.0;
	int size = width * height * 4;
	for (int i = 0; i < size; i += 4){
		grayVal = 0.299 * data[i] + 0.587 * data[i + 1] + 0.114 * data[i + 2];
		data[i] = data[i + 1] = data[i + 2] = (unsigned char)grayVal;
	}
	return true;
    //ClearToBlack();
    //return false;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	cout << "Quant Uniform" << endl;
	int size = width * height * 4;
	for (int i = 0; i < size; i += 4){
		data[i] /= 32;
		data[i] *= 32;
		data[i + 1] /= 32;
		data[i + 1] *= 32;
		data[i + 2] /= 64;
		data[i + 2] *= 64;
	}
	return true;
    //ClearToBlack();
    //return false;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	cout << "Quant Populosity" << endl;
	// data structure for histogram
	Color histo[256];
	int size = width * height * 4;

	// data structure for num of each color
	unsigned int RGB[32][32][32];
	// init data structure
	for (int i = 0; i < 32; ++i){
		for (int j = 0; j < 32; ++j){
			for (int k = 0; k < 32; ++k){
				RGB[i][j][k] = 0;
			}
		}
	}

	// quantize down to 32x32x32=32768 color values
	for (int i = 0; i < size; i += 4){
		data[i] /= 8;
		data[i + 1] /= 8;
		data[i + 2] /= 8;
	}

	// count num of each color 
	for (int i = 0; i < size; i += 4){
		//cout << (int)data[i];
		//cout << (int)data[i] << " " << (int)data[i + 1] << " " << (int)data[i + 2] << endl;
		++RGB[(int)data[i]][(int)data[i + 1]][(int)data[i + 2]];
	}

	// init search with first color
	Color currColor(data[0],data[1],data[2],RGB[(int)data[0]][(int)data[1]][(int)data[2]]);
	Color emptyColor(0, 0, 0, 0);
	for (int i = 0; i < 256; ++i){
		for (int j = 0; j < size; j += 4){
			if (RGB[(int)data[j]][(int)data[j + 1]][(int)data[j + 2]] > currColor.cnt){
				currColor.red = data[j];
				currColor.green = data[j + 1];
				currColor.blue = data[j + 2];
				currColor.cnt = RGB[(int)data[j]][(int)data[j + 1]][(int)data[j + 2]];
			}
		}
		histo[i] = currColor;
		RGB[currColor.red][currColor.green][currColor.blue] = 0;
		currColor = emptyColor;
	}
	
	// find closest color
	unsigned int currDist = 0;
	unsigned int currIdx = 0;
	unsigned int thisDist = 0;
	for (int i = 0; i < size; i += 4){
		currDist = (data[i] - histo[0].red)*(data[i] - histo[0].red) + (data[i + 1] - histo[0].green)*(data[i + 1] - histo[0].green) + (data[i + 2] - histo[0].blue)*(data[i + 2] - histo[0].blue);
		currIdx = 0;
		for (int j = 1; j < 256; ++j){
			thisDist = (data[i] - histo[j].red)*(data[i] - histo[j].red) + (data[i + 1] - histo[j].green)*(data[i + 1] - histo[j].green) + (data[i + 2] - histo[j].blue)*(data[i + 2] - histo[j].blue);
			if ( thisDist < currDist){
				currDist = thisDist;
				currIdx = j;
			}
		}
		data[i] = histo[currIdx].red;
		data[i + 1] = histo[currIdx].green;
		data[i + 2] = histo[currIdx].blue;
	}

	// re-map color values
	for (int i = 0; i < size; i += 4){
		data[i] *= 8;
		data[i + 1] *= 8;
		data[i + 2] *= 8;
	}

	return true;
    //ClearToBlack();
    //return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	cout << "Dither Threshold" << endl;
	int size = width * height * 4;
	float * imgData;
	imgData = new float[size];
	float grayVal = 0.0;

	// convert to grayscale
	for (int i = 0; i < size; i += 4){
		grayVal = 0.299 * data[i] + 0.587 * data[i + 1] + 0.114 * data[i + 2];
		data[i] = data[i + 1] = data[i + 2] = (unsigned char)grayVal;
	}

	// normalize
	for (int i = 0; i < size; i += 4){
		imgData[i] = (float)data[i] / 255;
		imgData[i+1] = (float)data[i+1] / 255;
		imgData[i+2] = (float)data[i+2] / 255;
	}

	// threshold
	for (int i = 0; i < size; i += 4){
		for (int j = i; j < i+3; ++j){
			(imgData[j] < 0.5) ? data[j] = 0 : data[j] = 255;
		}
	}

	delete [] imgData;
	return true;
    //ClearToBlack();
    //return false;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
	cout << "Dither Random" << endl;
	int size = width * height * 4;
	float * imgData;
	imgData = new float[size];
	float grayVal = 0.0;

	// convert to grayscale
	for (int i = 0; i < size; i += 4){
		grayVal = 0.299 * data[i] + 0.587 * data[i + 1] + 0.114 * data[i + 2];
		data[i] = (unsigned char)grayVal;
	}

	// normalize
	for (int i = 0; i < size; i += 4){
		imgData[i] = (float)data[i] / 255;
	}

	std::default_random_engine generator;
	double value;
	double mean = 0.0;
	double min = -0.2;
	double max = 0.2;
	double boundWidth = (-min + max) / 2;
	double stdDev = boundWidth / 3;
	std::normal_distribution<double> distribution(mean, stdDev);

	// add random noise from normal distribution
	for (int i = 0; i < size; i += 4){
		do{
			value = distribution(generator);
		} while (value <= min && value >= max);	
		imgData[i] += value;
	}

	// threshold
	for (int i = 0; i < size; i += 4){
		if (imgData[i] >= 0.5){
			data[i] = data[i + 1] = data[i + 2] = 255;
		}
		else{
			data[i] = data[i + 1] = data[i + 2] = 0;
		}

	}
	delete [] imgData;
	return true;
    //ClearToBlack();
    //return false;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
	cout << "Dither FS" << endl;
	int size = width * height * 4;
	float * imgData = new float[size];
	float grayVal = 0.0;

	// convert to grayscale
	for (int i = 0; i < size; i += 4){
		grayVal = 0.299 * data[i] + 0.587 * data[i + 1] + 0.114 * data[i + 2];
		data[i] = (unsigned char)grayVal;
	}

	// normalize
	for (int i = 0; i < size; i += 4){
		imgData[i] = (float)data[i] / 255;
	}

	int length = width * 4;
	int start = 0;
	int end = 0;
	float error = 0.0;
	float clusterMatrix[4] = { 0.4375, 0.0625, 0.3125, 0.1875 };
	for (int row = 0; row < height; ++row){
		if (row % 2 == 0){							// even
			start = row * length;					// start at begining of row
			end = start + length;					// end at end of row
			for (int i = start; i < end; i += 4){
				// determine error and threshold current pixel
				if (imgData[i] < 0.5){
					error = imgData[i];
					imgData[i] = 0.0;
				}
				else{
					error = imgData[i] - 1.0;
					imgData[i] = 1.0;
				}
				
				if (((i + 4) % length) != 0){  // checkout right bounds (i.e. not at end of row)
					imgData[i + 4] += error * clusterMatrix[0];
					if (row != (height-1)){				// check lower bounds (i.e. not on last row)
						imgData[i + length + 4] += error * clusterMatrix[1];
					}
				}
				if (row != (height - 1)){					// check lower bounds (i.e. not on last row)		
					imgData[i + length] += error * clusterMatrix[2];
					if ((i % length) != 0){					// check left bounds (i.e. not at begining of row)
						imgData[(i + length - 4)] += error * clusterMatrix[3];
					}
				}

			}
		}
		else{										// odd
			start = (row * length) + (length - 4);	// start at end of row
			end = start - length;					// end at begining of same row
			for (int i = start; i > end; i -= 4){
				// determine error and threshold current pixel
				if (imgData[i] < 0.5){
					error = imgData[i];
					imgData[i] = 0.0;
				}
				else{
					error = imgData[i] - 1.0;
					imgData[i] = 1.0;
				}

				if ((i % length) != 0){						// check left bounds (i.e. not at begining of row)
					imgData[i - 4] += error * clusterMatrix[0];
					if (row != (height - 1)){				// check lower bounds (i.e. not on last row)
						imgData[i + length - 4] += error * clusterMatrix[1];
					}
				}
				if (row != (height - 1)){					// check lower bounds (i.e. not on last row)		
					imgData[i + length] += error * clusterMatrix[2];
					if (((i + 4) % length) != 0){					// check right bounds (i.e. not at end of row)
						imgData[(i + length + 4)] += error * clusterMatrix[3];
					}
				}
			}
		}
	}

	// threshold, 0.5
	for (int i = 0; i < size; i += 4){
		if (imgData[i] >= 0.5){
			data[i] = data[i + 1] = data[i + 2] = 255;
		}
		else{
			data[i] = data[i + 1] = data[i + 2] = 0;
		}

	}
	delete [] imgData;
	return true;
    //ClearToBlack();
    //return false;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	cout << "Dither Bright" << endl;
	int size = width * height * 4;
	float * imgData = new float[size];
	float * sortedData = new float[size];
	float grayVal = 0.0;
	
	// convert to grayscale
	for (int i = 0; i < size; i += 4){
		grayVal = 0.299 * data[i] + 0.587 * data[i + 1] + 0.114 * data[i + 2];
		data[i] = (unsigned char)grayVal;
	}

	// normalize
	for (int i = 0; i < size; i += 4){
		imgData[i] = (float)data[i] / 255;
	}

	// copy to list for sorting
	for (int i = 0; i < size; i += 4){
		sortedData[i] = imgData[i];
	}

	// selection sort the colors;
	int currMin = 0;
	float tmpColor = 0.0;
	for (int i = 0; i < size-4; i+=4){
		currMin = i;
		for (int j = i + 4; j < size; j+=4){
			if (sortedData[j] < sortedData[currMin]){
				currMin = j;
			}
		}
		if (currMin != i){
			tmpColor = sortedData[i];
			sortedData[i] = sortedData[currMin];
			sortedData[currMin] = tmpColor;
		}
	}

	// threshold
	float threshold = sortedData[size / 2];
	//cout << "threshold: " << threshold << endl;
	for (int i = 0; i < size; i += 4){
		if (imgData[i] >= threshold){
			data[i] = data[i + 1] = data[i + 2] = 255;
		}
		else{
			data[i] = data[i + 1] = data[i + 2] = 0;
		}

	}

	delete [] imgData;
	delete [] sortedData;
	return true;

    //ClearToBlack();
    //return false;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
	cout << "Dither Cluster" << endl;
	int size = width * height * 4;
	float * imgData = new float[size];
	float grayVal = 0.0;
	float clusterMatrix[16] = { 0.7500, 0.3750, 0.6250, 0.2500,
								0.0625, 1.0000, 0.8750, 0.4375,
								0.5000, 0.8125, 0.9375, 0.1250,
								0.1875, 0.5625, 0.3125, 0.6875 };

	// convert to grayscale
	for (int i = 0; i < size; i += 4){
		grayVal = 0.299 * data[i] + 0.587 * data[i + 1] + 0.114 * data[i + 2];
		data[i] = (unsigned char)grayVal;
	}

	// normalize
	for (int i = 0; i < size; i += 4){
		imgData[i] = (float)data[i] / 255;
	}

	int idx = 0;
	for (int i = 0; i < size - 16; i += 16){
		// at the begining of each row(except first), skip 3 rows, they are already processed
		if ((i != 0) && (i % (width * 4) == 0)){
			// skip 3 rows
			i += ((width * 4) * 3);
		}
		
		// set start of cluster
		idx = i;
		// walk through 16 cluster thresholds
		for (int j = 0; (j < 16) && (idx < size) ; ++j){
			// threshold
			if (imgData[idx] >= clusterMatrix[j]){
				data[idx] = data[idx + 1] = data[idx + 2] = 255;
			}
			else{
				data[idx] = data[idx + 1] = data[idx + 2] = 0;
			}
			// skip to next row of data at the end of each row of cluster matrix
			if ((j % 4) == 3){
				if (j != 15){
					idx += (width * 4) - 12;
				}
			}
			else{
				idx += 4;
			}
		}

	}
	
	delete [] imgData;
	return true;

    //ClearToBlack();
    //return false;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
	cout << "Dither Color" << endl;
	int size = width * height * 4;
	double * imgData = new double[size];
	//float * colorData = new float[size];
	float colorTable8[8];
	colorTable8[0] = 0.0; colorTable8[1] = 36.0 / 255; colorTable8[2] = 73.0 / 255; colorTable8[3] = 109.0 / 255;
	colorTable8[4] = 146.0 / 255; colorTable8[5] = 182.0 / 255; colorTable8[6] = 219.0 / 255; colorTable8[7] = 1.0;
	float colorTable4[4];
	colorTable4[0] = 0.0; colorTable4[1] = 85.0 / 255; colorTable4[2] = 170.0 / 255; colorTable4[3] = 1.0;
	//float grayVal = 0.0;

//	for (int i = 0; i < size; i += 4){
//		data[i] /= 32;
//		data[i + 1] /= 32;
//		data[i + 2] /= 64;
//	}

	// normalize
	for (int i = 0; i < size; i += 4){
		imgData[i] = (double)data[i] / 255;
		imgData[i+1] = (double)data[i+1] / 255;
		imgData[i+2] = (double)data[i+2] / 255;
	}

	int length = width * 4;
	int start = 0;
	int end = 0;
	ColorFT error;
	float thisDist = 0.0;
	float currDist = 0.0;
	float clusterMatrix[4] = { 0.4375, 0.0625, 0.3125, 0.1875 };
	for (int row = 0; row < height; ++row){
		if (row % 2 == 0){							// even
			start = row * length;					// start at begining of row
			end = start + length;					// end at end of row
			//cout << "even " << start << ' ' << end <<  endl;
			for (int i = start; i < end; i += 4){
				// determine error 
				thisDist = 0.0;
				currDist = 999.0;
				for (int red = 0; red < 8; ++red){
					for (int green = 0; green < 8; ++green){
						for (int blue = 0; blue < 4; ++blue){
							thisDist = sqrt((imgData[i] - colorTable8[red])*(imgData[i] - colorTable8[red]) + (imgData[i + 1] - colorTable8[green])*(imgData[i + 1] - colorTable8[green]) + (imgData[i + 2] - colorTable4[blue])*(imgData[i + 2] - colorTable4[blue]));
							if (thisDist < currDist){
								currDist = thisDist;
								error.red = colorTable8[red];
								error.green = colorTable8[green];
								error.blue = colorTable4[blue];
							}
						}
					}
				}
				// set current pixel's new color value
				imgData[i] = error.red;
				imgData[i + 1] = error.green;
				imgData[i + 2] = error.blue;

				if (((i + 4) % length) != 0){  // checkout right bounds (i.e. not at end of row)
					imgData[i + 4] = min(1.0, max(0.0, imgData[i+4] + error.red * clusterMatrix[0]));
					imgData[i + 5] = min(1.0, max(0.0, imgData[i+5] + error.green * clusterMatrix[0]));
					imgData[i + 6] = min(1.0, max(0.0, imgData[i+6] + error.blue * clusterMatrix[0]));
					if (row != (height - 1)){				// check lower bounds (i.e. not on last row)
						imgData[i + length + 4] = min(1.0, max(0.0, imgData[i + length + 4] + error.red * clusterMatrix[1]));
						imgData[i + length + 5] = min(1.0, max(0.0,imgData[i + length + 5] + error.green * clusterMatrix[1]));
						imgData[i + length + 6] = min(1.0, max(0.0, imgData[i + length + 6] + error.blue * clusterMatrix[1]));
					}
				}
				if (row != (height - 1)){					// check lower bounds (i.e. not on last row)		
					imgData[i + length] = min(1.0, max(0.0, imgData[i + length] + error.red * clusterMatrix[2]));
					imgData[i + length + 1] = min(1.0, max(0.0,imgData[i + length + 1] + error.green * clusterMatrix[2]));
					imgData[i + length + 2] = min(1.0, max(0.0, imgData[i + length + 2] + error.blue * clusterMatrix[2]));
					if ((i % length) != 0){					// check left bounds (i.e. not at begining of row)
						imgData[(i + length - 4)] = min(1.0, max(0.0, imgData[(i + length - 4)] + error.red* clusterMatrix[3]));
						imgData[(i + length - 3)] = min(1.0, max(0.0, imgData[(i + length - 3)] + error.green * clusterMatrix[3]));
						imgData[(i + length - 2)] = min(1.0, max(0.0, imgData[(i + length - 2)] + error.blue * clusterMatrix[3]));
					}
				}

			}
		}
		else{										// odd
			start = (row * length) + (length - 4);	// start at end of row
			end = start - length;					// end at begining of same row
			//cout << "odd " << start << ' ' << end << endl;
			for (int i = start; i > end; i -= 4){
				// determine error 
				thisDist = 0.0;
				currDist = 999.0;
				for (int red = 0; red < 8; ++red){
					for (int green = 0; green < 8; ++green){
						for (int blue = 0; blue < 4; ++blue){
							thisDist = sqrt((imgData[i] - colorTable8[red])*(imgData[i] - colorTable8[red]) + (imgData[i + 1] - colorTable8[green])*(imgData[i + 1] - colorTable8[green]) + (imgData[i + 2] - colorTable4[blue])*(imgData[i + 2] - colorTable4[blue]));
							if (thisDist < currDist){
								currDist = thisDist;
								error.red = colorTable8[red];
								error.green = colorTable8[green];
								error.blue = colorTable4[blue];
							}
						}
					}
				}
				// set current pixel's new color value
				imgData[i] = error.red;
				imgData[i + 1] = error.green;
				imgData[i + 2] = error.blue;

				if ((i % length) != 0){						// check left bounds (i.e. not at begining of row)
					imgData[i - 4] = min(1.0, max(0.0, imgData[i - 4] + error.red * clusterMatrix[0]));
					imgData[i - 3] = min(1.0, max(0.0, imgData[i - 3] + error.green * clusterMatrix[0]));
					imgData[i - 2] = min(1.0, max(0.0, imgData[i - 2] + error.blue * clusterMatrix[0]));
					if (row != (height - 1)){				// check lower bounds (i.e. not on last row)
						imgData[i + length - 4] = min(1.0, max(0.0, imgData[i + length - 4] + error.red * clusterMatrix[1]));
						imgData[i + length - 3] = min(1.0, max(0.0, imgData[i + length - 3] + error.green * clusterMatrix[1]));
						imgData[i + length - 2] = min(1.0, max(0.0, imgData[i + length - 2] + error.blue * clusterMatrix[1]));
					}
				}
				if (row != (height - 1)){					// check lower bounds (i.e. not on last row)		
					imgData[i + length] = min(1.0, max(0.0, imgData[i + length] + error.red * clusterMatrix[2]));
					imgData[i + length + 1] = min(1.0, max(0.0, imgData[i + length + 1] + error.green * clusterMatrix[2]));
					imgData[i + length + 2] = min(1.0, max(0.0, imgData[i + length + 2] + error.blue * clusterMatrix[2]));
					if (((i + 4) % length) != 0){					// check right bounds (i.e. not at end of row)
						imgData[(i + length + 4)] = min(1.0, max(0.0, imgData[(i + length + 4)] + error.red * clusterMatrix[3]));
						imgData[(i + length + 5)] = min(1.0, max(0.0, imgData[(i + length + 5)] + error.green * clusterMatrix[3]));
						imgData[(i + length + 6)] = min(1.0, max(0.0, imgData[(i + length + 6)] + error.blue * clusterMatrix[3]));
					}
				}
			}
		}
	}

	// remap to full color value
	for (int i = 0; i < size; i += 4){
		data[i] = imgData[i] * 255;
		data[i + 1] = imgData[i+1] * 255;
		data[i + 2] = imgData[i+2] * 255;
	}

	delete [] imgData;
	return true;
    //ClearToBlack();
    //return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
	cout << "Comp Over" << endl;
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }


	int size = width * height * 4;
	//float alphaF = 0.05;
	//float alphaG = 1.0-alphaF;
	float F = 1.0;
	float G;

	// add foreground and background
	for (int i = 0; i < size; i += 4){
		G = 1 - (float)data[i + 3] / 255;
		data[i] = data[i] * F + pImage->data[i] * G;
		data[i+1] = data[i+1] * F + pImage->data[i+1] * G;
		data[i+2] = data[i+2] * F + pImage->data[i+2] * G;
		data[i + 3] = data[i+3] * F  + pImage->data[i+3] * G;
	}

	return true;
    //ClearToBlack();
    //return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
	cout << "Comp In" << endl;
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

	int size = width * height * 4;
	float F = 1.0;
	float G = 0;

	// add foreground and background
	for (int i = 0; i < size; i += 4){
		F = (float)pImage->data[i + 3] / 255;
		data[i] = data[i] * F;
		data[i + 1] = data[i + 1] * F;
		data[i + 2] = data[i + 2] * F;
		data[i + 3] = data[i + 3] * F;
	}

	return true;
    //ClearToBlack();
    //return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
	cout << "Comp Out" << endl;
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

	int size = width * height * 4;
	float F = 1.0;
	float G = 0;

	// add foreground and background
	for (int i = 0; i < size; i += 4){
		F = 1 - (float)pImage->data[i + 3] / 255;
		data[i] = data[i] * F;
		data[i + 1] = data[i + 1] * F;
		data[i + 2] = data[i + 2] * F;
		data[i + 3] = data[i + 3] * F;
	}

	return true;

    //ClearToBlack();
    //return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
	cout << "Comp Atop" << endl;
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

	int size = width * height * 4;
	float F = 1.0;
	float G = 0;
	
	// add foreground and background
	for (int i = 0; i < size; i += 4){
		F = (float)pImage->data[i + 3] / 255;
		G = 1 - (float)data[i + 3] / 255;
		data[i] = data[i] * F + pImage->data[i] * G;
		data[i + 1] = data[i + 1] * F + pImage->data[i + 1] * G;
		data[i + 2] = data[i + 2] * F + pImage->data[i + 2] * G;
		data[i + 3] = data[i + 3] * F + pImage->data[i + 3] * G;
	}

	return true;
    //ClearToBlack();
    //return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
	cout << "Comp Xor" << endl;
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }


	int size = width * height * 4;
	float F;
	float G;

	// add foreground and background
	for (int i = 0; i < size; i += 4){
		F = 1.0 - (float)pImage->data[i + 3]/255;
		G = 1.0 - (float)data[i + 3]/255;
		data[i] = data[i] * F + pImage->data[i] * G;
		data[i + 1] = data[i + 1] * F + pImage->data[i + 1] * G;
		data[i + 2] = data[i + 2] * F + pImage->data[i + 2] * G;
		data[i + 3] = data[i + 3] * F + pImage->data[i + 3] * G;
	}

    //ClearToBlack();
    //return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	cout << "Filter Box" << endl;
	// init 2D array for Gaussian Filter
	//float gaussWeight = 1.0 / 256;
	float boxFilter[5][5] =
	{ { 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25 },
	{ 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25 },
	{ 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25 },
	{ 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25 },
	{ 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25, 1.0 / 25 } };

	// after filter image will be smaller
	int newHeight = height - 4;
	int newWidth = width - 4;

	// init 2D array of colors data struct
	Color** imgData = new Color*[height];
	for (int i = 0; i < height; ++i){
		imgData[i] = new Color[width];
	}

	// init 2D array of colors data struct
	Color** filtData = new Color*[newHeight];
	for (int i = 0; i < newHeight; ++i){
		filtData[i] = new Color[newWidth];
	}
	for (int row = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col){
			filtData[row][col].red = 0;
			filtData[row][col].green = 0;
			filtData[row][col].blue = 0;
		}
	}

	//	load raw pixel data into 2D array of colors
	int pos = 0;
	for (int row = 0; row < height; ++row){
		for (int col = 0; col < width; ++col, pos += 4){
			imgData[row][col].red = data[pos];
			imgData[row][col].green = data[pos + 1];
			imgData[row][col].blue = data[pos + 2];
		}
	}

	for (int row = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col){
			for (int i = 0; i < 5; ++i){
				for (int j = 0; j < 5; ++j){
					// check pixels are in bound
					//if (row + i >= 0 && col + j >= 0 && row + i < newHeight && col + j < newWidth){
					filtData[row][col].red += imgData[row + i][col + j].red * boxFilter[i][j];
					filtData[row][col].green += imgData[row + i][col + j].green * boxFilter[i][j];
					filtData[row][col].blue += imgData[row + i][col + j].blue * boxFilter[i][j];
					//}
				}
			}
		}
	}


	// put 2D array of colors back into 1D data array
	for (int row = 0, pos = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col, pos += 4){
			data[pos] = filtData[row][col].red;
			data[pos + 1] = filtData[row][col].green;
			data[pos + 2] = filtData[row][col].blue;
		}
	}

	// clean up memory
	for (int i = 0; i < height; ++i){
		delete[] imgData[i];
	}
	delete[] imgData;

	for (int i = 0; i < newHeight; ++i){
		delete[] filtData[i];
	}
	delete[] filtData;

	// set new dims of output image
	width = newWidth;
	height = newHeight;

	return true;
    //ClearToBlack();
    //return false;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	cout << "Filter Bartlett" << endl;
	// init 2D array for Gaussian Filter
	float bartFilter[5][5] =
	{ { 1.0 / 81, 2.0 / 81, 3.0 / 81, 2.0 / 81, 1.0 / 81 },
	{ 2.0 / 81, 4.0 / 81, 6.0 / 81, 4.0 / 81, 2.0 / 81 },
	{ 3.0 / 81, 6.0 / 81, 9.0 / 81, 6.0 / 81, 3.0 / 81 },
	{ 2.0 / 81, 4.0 / 81, 6.0 / 81, 4.0 / 81, 2.0 / 81 },
	{ 1.0 / 81, 2.0 / 81, 3.0 / 81, 2.0 / 81, 1.0 / 81 } };

	// after filter image will be smaller
	int newHeight = height - 4;
	int newWidth = width - 4;

	// init 2D array of colors data struct
	Color** imgData = new Color*[height];
	for (int i = 0; i < height; ++i){
		imgData[i] = new Color[width];
	}

	// init 2D array of colors data struct
	Color** filtData = new Color*[newHeight];
	for (int i = 0; i < newHeight; ++i){
		filtData[i] = new Color[newWidth];
	}
	for (int row = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col){
			filtData[row][col].red = 0;
			filtData[row][col].green = 0;
			filtData[row][col].blue = 0;
		}
	}

	//	load raw pixel data into 2D array of colors
	int pos = 0;
	for (int row = 0; row < height; ++row){
		for (int col = 0; col < width; ++col, pos += 4){
			imgData[row][col].red = data[pos];
			imgData[row][col].green = data[pos + 1];
			imgData[row][col].blue = data[pos + 2];
		}
	}

	for (int row = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col){
			for (int i = 0; i < 5; ++i){
				for (int j = 0; j < 5; ++j){
					// check pixels are in bound
					//if (row + i >= 0 && col + j >= 0 && row + i < newHeight && col + j < newWidth){
					filtData[row][col].red += imgData[row + i][col + j].red * bartFilter[i][j];
					filtData[row][col].green += imgData[row + i][col + j].green * bartFilter[i][j];
					filtData[row][col].blue += imgData[row + i][col + j].blue * bartFilter[i][j];
					//}
				}
			}
		}
	}


	// put 2D array of colors back into 1D data array
	for (int row = 0, pos = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col, pos += 4){
			data[pos] = filtData[row][col].red;
			data[pos + 1] = filtData[row][col].green;
			data[pos + 2] = filtData[row][col].blue;
		}
	}

	// clean up memory
	for (int i = 0; i < height; ++i){
		delete[] imgData[i];
	}
	delete[] imgData;

	for (int i = 0; i < newHeight; ++i){
		delete[] filtData[i];
	}
	delete[] filtData;

	// set new dims of output image
	width = newWidth;
	height = newHeight;

	return true;
    //ClearToBlack();
    //return false;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	cout << "Filter Gaussian" << endl;
	// init 2D array for Gaussian Filter
	float gaussWeight = 1.0 / 256;
	float gaussFilter[5][5] = 
	{ { 1.0/256, 4.0/256, 6.0/256, 4.0/256, 1.0/256 },
	{ 4.0/256, 16.0/256, 24.0/256, 16.0/256, 4.0/256 },
	{ 6.0/256, 24.0/256, 36.0/256, 24.0/256, 6.0/256 },
	{ 4.0/256, 16.0/256, 24.0/256, 16.0/256, 4.0/256 },
	{ 1.0/256, 4.0/256, 6.0/256, 4.0/256, 1.0/256 } };

	// after filter image will be smaller
	int newHeight = height - 4;
	int newWidth = width - 4;
		
	// init 2D array of colors data struct
	Color** imgData = new Color*[height];
	for (int i = 0; i < height; ++i){
		imgData[i] = new Color[width];
	}

	// init 2D array of colors data struct
	Color** filtData = new Color*[newHeight];
	for (int i = 0; i < newHeight; ++i){
		filtData[i] = new Color[newWidth];
	}
	for (int row = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col){
			filtData[row][col].red = 0;
			filtData[row][col].green = 0;
			filtData[row][col].blue = 0;
		}
	}

	//	load raw pixel data into 2D array of colors
	int pos = 0;
	for (int row = 0; row < height; ++row){
		for (int col = 0; col < width; ++col, pos+=4){
			imgData[row][col].red = data[pos];
			imgData[row][col].green = data[pos + 1];
			imgData[row][col].blue = data[pos + 2];
		}
	}

	for (int row = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col){
			for (int i = 0; i < 5; ++i){
				for (int j = 0; j < 5; ++j){
					// check pixels are in bound
					//if (row + i >= 0 && col + j >= 0 && row + i < newHeight && col + j < newWidth){
						filtData[row][col].red += imgData[row + i][col + j].red * gaussFilter[i][j];
						filtData[row][col].green += imgData[row + i][col + j].green * gaussFilter[i][j];
						filtData[row][col].blue += imgData[row + i][col + j].blue * gaussFilter[i][j];
					//}
				}
			}
		}
	}


	// put 2D array of colors back into 1D data array
	for (int row = 0, pos = 0; row < newHeight; ++row){
		for (int col = 0; col < newWidth; ++col, pos += 4){
			data[pos] = filtData[row][col].red;
			data[pos + 1] = filtData[row][col].green;
			data[pos + 2] = filtData[row][col].blue;
		}
	}

	// clean up memory
	for (int i = 0; i < height; ++i){
		delete [] imgData[i];
	}
	delete [] imgData;

	for (int i = 0; i < newHeight; ++i){
		delete[] filtData[i];
	}
	delete [] filtData;

	// set new dims of output image
	width = newWidth;
	height = newHeight;

	return true;
    //ClearToBlack();
    //return false;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	int size = height/2 * width * 4;
	for (int i = 0; i < size; i += 4){
		data[i] = 0;
		data[i + 1] = 0;
		data[i + 2] = 0;
		data[i + 3] = 0;

	}
	return true;
    //ClearToBlack();
    //return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

