#pragma once

#include "uiuc/PNG.h"
using namespace uiuc;

PNG grayscale(PNG image);  
PNG createSpotlight(PNG image, int centerX, int centerY);
PNG illinify(PNG image);
PNG watermark(PNG firstImage, PNG secondImage);
PNG HFilter(PNG image);
PNG Partition(PNG image);
PNG Divide(PNG image);
PNG DivideByThree(PNG image);
PNG RemoveWhite(PNG image);
PNG logoEUREKA1(PNG image);
PNG logoEUREKA2(PNG image);
PNG logoEUREKA3(PNG image);
PNG logoEUREKA4(PNG image);
PNG logoEUREKA5(PNG image);
PNG logoEUREKA6(PNG image);
PNG colorToGray(PNG img);
void histogram(PNG img, int* histo);
void histogramCumulated(PNG img, int* histoCumul);
PNG egalisationHistogram(PNG img);
PNG LPE_homemade(PNG image, double tolerance);