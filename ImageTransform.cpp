#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "uiuc/PNG.h"
#include "uiuc/HSLAPixel.h"
#include "ImageTransform.h"
#include "Point.hpp"
#include "Pt.hpp"

/* ******************
(Begin multi-line comment...)

Write your name and email push_backress in the comment space here:

Name:
Email:

(...end multi-line comment.)
******************** */

using uiuc::PNG;
using uiuc::HSLAPixel;

/**
 * Returns an image that has been transformed to grayscale.
 *
 * The saturation of every pixel is set to 0, removing any color.
 *
 * @return The grayscale image.
 */
PNG grayscale(PNG image) {
  /// This function is already written for you so you can see how to
  /// interact with our PNG class.
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);

      // `pixel` is a reference to the memory stored inside of the PNG `image`,
      // which means you're changing the image directly. No need to `set`
      // the pixel since you're directly changing the memory of the image.
      pixel.s = 0;
    }
  }

  return image;
}



/**
 * Returns an image with a spotlight centered at (`centerX`, `centerY`).
 *
 * A spotlight adjusts the luminance of a pixel based on the distance the pixel
 * is away from the center by decreasing the luminance by 0.5% per 1 pixel euclidean
 * distance away from the center.
 *
 * For example, a pixel 3 pixels above and 4 pixels to the right of the center
 * is a total of `sqrt((3 * 3) + (4 * 4)) = sqrt(25) = 5` pixels away and
 * its luminance is decreased by 2.5% (0.975x its original value).  At a
 * distance over 160 pixels away, the luminance will always decreased by 80%.
 * 
 * The modified PNG is then returned.
 *
 * @param image A PNG object which holds the image data to be modified.
 * @param centerX The center x coordinate of the crosshair which is to be drawn.
 * @param centerY The center y coordinate of the crosshair which is to be drawn.
 *
 * @return The image with a spotlight.
 */
PNG createSpotlight(PNG image, int centerX, int centerY) {
int xx=0, yy=0,z;	
  for (unsigned x = 0; x < image.width(); x++) {
        for (unsigned y = 0; y < image.height(); y++) {
        HSLAPixel & pixel = image.getPixel(x, y);
		xx=x-centerX;
		yy=y-centerY;
		z=sqrt((xx*xx)+(yy*yy));
		if(z<160)
       	        pixel.l=pixel.l-.005*z*pixel.l;
		else
                pixel.l = pixel.l - 0.8*pixel.l;

        }
    }
  return image;
  
}
 

/**
 * Returns a image transformed to Illini colors.
 *
 * The hue of every pixel is set to the a hue value of either orange or
 * blue, based on if the pixel's hue value is closer to orange than blue.
 *
 * @param image A PNG object which holds the image data to be modified.
 *
 * @return The illinify'd image.
**/
PNG illinify(PNG image) {
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);

      if (pixel.h>=101&&pixel.h<281)
      pixel.h = 216 ;
      else
	pixel.h =11 ;
	
    }
  }

  return image;
}
 

/**
* Returns an immge that has been watermarked by another image.
*
* The luminance of every pixel of the second image is checked, if that
* pixel's luminance is 1 (100%), then the pixel at the same location on
* the first image has its luminance increased by 0.2.
*
* @param firstImage  The first of the two PNGs, which is the base image.
* @param secondImage The second of the two PNGs, which acts as the stencil.
*
* @return The watermarked image.
*/
PNG watermark(PNG firstImage, PNG secondImage) {
  for (unsigned x = 0; x < firstImage.width(); x++) {
    for (unsigned y = 0; y < firstImage.height(); y++) {
      HSLAPixel & pixel1 = firstImage.getPixel(x, y);
      HSLAPixel & pixel2 = secondImage.getPixel(x, y);
	if(pixel2.l==1&&pixel1.l+.2<1)
	pixel1.l=pixel1.l+.2;	
	
	}
}
  return firstImage;
}

PNG HFilter(PNG image) {
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);
    }
  }

  return image;
}

PNG Divide(PNG image) {
  uiuc::PNG imgout(image.width()/2,image.height()/2);

  for (unsigned x = 0; x < imgout.width(); x++) {
    for (unsigned y = 0; y < imgout.height(); y++) {
      HSLAPixel & pixel = imgout.getPixel(x, y);

      HSLAPixel & pixel1 = image.getPixel(x*2, y*2);
      HSLAPixel & pixel2 = image.getPixel(x*2+1, y*2);
      HSLAPixel & pixel3 = image.getPixel(x*2, y*2+1);
      HSLAPixel & pixel4 = image.getPixel(x*2+1, y*2+1);

      pixel.h = (pixel1.h + pixel2.h + pixel3.h + pixel4.h)/4; //over 360
      pixel.s = (pixel1.s + pixel2.s + pixel3.s + pixel4.s)/4; //over 1
      pixel.l = (pixel1.l + pixel2.l + pixel3.l + pixel4.l)/4; //over 1  
      pixel.a = (pixel1.a + pixel2.a + pixel3.a + pixel4.a)/4;
    }
  }

  return imgout;
}

PNG DivideByThree(PNG image) {
  uiuc::PNG imgout(image.width()/3,image.height()/3);

  for (unsigned x = 0; x < imgout.width(); x++) {
    for (unsigned y = 0; y < imgout.height(); y++) {
      HSLAPixel & pixel = imgout.getPixel(x, y);

      HSLAPixel & pixel1 = image.getPixel(x*3, y*3);
      HSLAPixel & pixel2 = image.getPixel(x*3+1, y*3);
      HSLAPixel & pixel3 = image.getPixel(x*3+2, y*3);

      HSLAPixel & pixel4 = image.getPixel(x*3, y*3+1);
      HSLAPixel & pixel5 = image.getPixel(x*3+1, y*3+1);
      HSLAPixel & pixel6 = image.getPixel(x*3+2, y*3+1);

      HSLAPixel & pixel7 = image.getPixel(x*3, y*3+2);
      HSLAPixel & pixel8 = image.getPixel(x*3+1, y*3+2);
      HSLAPixel & pixel9 = image.getPixel(x*3+2, y*3+2);

      pixel.h = (pixel1.h + pixel2.h + pixel3.h + pixel4.h + pixel5.h + pixel6.h + pixel7.h + pixel8.h + pixel9.h)/9; //over 360
      pixel.s = (pixel1.s + pixel2.s + pixel3.s + pixel4.s + pixel5.s + pixel6.s + pixel7.s + pixel8.s + pixel9.s)/9; //over 1
      pixel.l = (pixel1.l + pixel2.l + pixel3.l + pixel4.l + pixel5.l + pixel6.l + pixel7.l + pixel8.l + pixel9.l)/9; //over 1  
      pixel.a = (pixel1.a + pixel2.a + pixel3.a + pixel4.a + pixel5.a + pixel6.a + pixel7.a + pixel8.a + pixel9.a)/9; //over 1
    }
  }

  return imgout;
}

PNG RemoveWhite(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);

      if(pixel.s <= 15 && pixel.l >= 0.70){
        pixel.a = 0;
      }
    }
  }

  return image;
}

PNG logoEUREKA1(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);

      if(pixel.h <= 30 || pixel.l <= 0.45){
        pixel.l = 1;
      }
    }
  }

  return image;
}

PNG logoEUREKA2(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);
      pixel.l = 0.50;
      pixel.s = 0;
    }
  }

  return image;
}

PNG logoEUREKA3(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);
      pixel.l = 1;
      pixel.s = 0;
    }
  }

  return image;
}

PNG logoEUREKA4(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);
      pixel.l = 0;
      pixel.s = 0;
    }
  }

  return image;
}

PNG logoEUREKA5(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);
      pixel.h = 306;
      pixel.s = 0.81;
      pixel.l = 0.59;
    }
  }

  return image;
}

PNG Partition(PNG image) {
  PNG Sol_key, Fa_key, Curly_bracket;
  Sol_key.readFromFile("SolKey.png");
  Fa_key.readFromFile("FaKey.png");
  Curly_bracket.readFromFile("CurlyBracket.png");

  //White image
  for(unsigned int i = 0; i < image.height(); i++){
    for(unsigned int j = 0; j < image.width(); j++){
      HSLAPixel & pixel = image.getPixel(j,i);
      pixel.h = 0; //over 360
      pixel.s = 0; //over 1
      pixel.l = 1; //over 1  
      pixel.a = 1;
    }
  }

  unsigned margin_LEFT = 200;
  unsigned margin_RIGHT = 150;
  unsigned margin_UP = 600;
  unsigned margin_BETWEEN_LINES = 150;
  unsigned margin_DOWN = 250;
  unsigned interval = 20;
  unsigned margin_SPECIAL = 500;

  Point l1;
  Point l2(0,interval);
  Point l3(0,interval*2);
  Point l4(0,interval*3);
  Point l5(0,interval*4);

  Point l6(0,interval*10);
  Point l7(0,interval*11);
  Point l8(0,interval*12);
  Point l9(0,interval*13);
  Point l10(0,interval*14);

  //Draw black 5 lines
  for (unsigned x = margin_LEFT; x < image.width() - margin_RIGHT; x++) {
    for (unsigned y = margin_UP; y < image.height() - margin_DOWN; y++) {
      if( 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l1.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l2.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l3.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l4.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l5.getY() ||
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l6.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l7.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l8.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l9.getY() || 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) == l10.getY()
      ) {
        HSLAPixel & pixel = image.getPixel(x, y);
        pixel.h = 0; //over 360
        pixel.s = 0; //over 1
        pixel.l = 0; //over 1
        pixel.a = 1;
      }
      else{
        HSLAPixel & pixel = image.getPixel(x, y);
        pixel.h = 0; //over 360
        pixel.s = 0; //over 1
        pixel.l = 1; //over 1  
        pixel.a = 1;
      }
    }
  }
  

  //Draw vertical lines
  for (unsigned x = margin_LEFT; x < image.width() - margin_RIGHT; x++) {
    for (unsigned y = margin_UP; y < image.height() - margin_DOWN; y++) {
      if( 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) <= l10.getY() && 
        (
          (x-margin_LEFT == 0) || //First vertical lines
          ( (x-margin_LEFT-margin_SPECIAL)%((image.width()-margin_RIGHT-margin_LEFT-margin_SPECIAL)/4) == 0 && x >= (margin_LEFT + margin_SPECIAL) ) || //Others vertical lines
          ( x >= margin_LEFT-Curly_bracket.width() && x < margin_LEFT )//&& Curly_bracket.getPixel(x-margin_LEFT+Curly_bracket.width(),(y-margin_UP)%(margin_BETWEEN_LINES+interval*15)).a != 0 )
        )
      ) {
        HSLAPixel & pixel = image.getPixel(x, y);
        pixel.h = 0; //over 360
        pixel.s = 0; //over 1
        pixel.l = 0; //over 1
        pixel.a = 1; //over 1
      }
    }
  }
  
  //Draw Curly Brackets
  for (unsigned x = 0; x < margin_LEFT; x++) {
    for (unsigned y = margin_UP; y < image.height() - margin_DOWN; y++) {
      if( 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) <= l10.getY()+24 && 
        ( x >= margin_LEFT-Curly_bracket.width() && x < margin_LEFT && Curly_bracket.getPixel(x-margin_LEFT+Curly_bracket.width(),(y-margin_UP)%(margin_BETWEEN_LINES+interval*15)).a != 0 )
      ){
        HSLAPixel & pixel = image.getPixel(x, y-13);
        pixel.h = 0; //over 360
        pixel.s = 0; //over 1
        pixel.l = 0; //over 1
        pixel.a = 1; //over 1
      }
    }
  }

  //Draw Sol keys
  for (unsigned x = margin_LEFT; x < margin_LEFT + Sol_key.width(); x++) {
    for (unsigned y = margin_UP-50; y < image.height() - margin_DOWN; y++) {
      if( 
        (y-margin_UP-50)%(margin_BETWEEN_LINES+interval*15) < Sol_key.height() && 
        ( x <= margin_LEFT+Sol_key.width() && Sol_key.getPixel(x-margin_LEFT,(y-margin_UP-50)%(margin_BETWEEN_LINES+interval*15)).a != 0 )
      ){
        HSLAPixel & pixel = image.getPixel(x, y-90);
        pixel.h = 0; //over 360
        pixel.s = 0; //over 1
        pixel.l = 0; //over 1
        pixel.a = 1; //over 1
      }
    }
  }

  //Draw Fa keys
  for (unsigned x = margin_LEFT; x < margin_LEFT + Fa_key.width(); x++) {
    for (unsigned y = margin_UP; y < image.height() - margin_DOWN; y++) {
      if( 
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) >= l6.getY() &&
        (y-margin_UP)%(margin_BETWEEN_LINES+interval*15) < l6.getY() + Fa_key.height() && 
        (Fa_key.getPixel(x-margin_LEFT,(y-margin_UP-l6.getY())%(margin_BETWEEN_LINES+interval*15)).a != 0 )
      ){
        HSLAPixel & pixel = image.getPixel(x, y-10);
        pixel.h = 0; //over 360
        pixel.s = 0; //over 1
        pixel.l = 0; //over 1
        pixel.a = 1; //over 1
      }
    }
  }

  return image;
}

PNG logoEUREKA6(PNG image){
  for (unsigned x = 0; x < image.width(); x++) {
    for (unsigned y = 0; y < image.height(); y++) {
      HSLAPixel & pixel = image.getPixel(x, y);
      pixel.h = y%360;

      if(pixel.l >= 0.45){
        pixel.l = 0.60;
      }
    }
  }

  return image;
}

PNG colorToGray(PNG img){
  uiuc::PNG imgout(img.width(),img.height());

  for (unsigned x = 0; x < img.width(); x++) {
    for (unsigned y = 0; y < img.height(); y++) {
      HSLAPixel & pixel = img.getPixel(x, y);
      HSLAPixel & pixel_ = imgout.getPixel(x, y);

      pixel_.h = pixel.h;
      pixel.s = 0;
      pixel_.l = pixel.l;
      pixel_.a = pixel.a;
    }
  }

  return imgout;
}

void histogram(PNG img, int* histo){
  int gray_level;

  for (unsigned x = 0; x < img.width(); x++) {
    for (unsigned y = 0; y < img.height(); y++) {
      HSLAPixel & pixel = img.getPixel(x, y);
      gray_level = (double) (pixel.l * 255);
      histo[gray_level]++;
    }
  }
}

void histogramCumulated(PNG img, int* histoCumul){//size = 256 for histoCumul
  int histo[256] = {};
  histogram(img, histo);

  for(int k = 0; k < 256; k++){
    for(int i = 0; i <= k; i++){
      histoCumul[k] += histo[i];
    }
  }
}

PNG egalisationHistogram(PNG img){
  uiuc::PNG imgout = colorToGray(img);
  int C[256] = {};
  int N = img.height() * img.width();
  histogramCumulated(img,C);
  int gray_level;
  
  for (unsigned x = 0; x < img.width(); x++) {
    for (unsigned y = 0; y < img.height(); y++) {
      HSLAPixel & pixel = img.getPixel(x, y);
      HSLAPixel & pixel_ = imgout.getPixel(x, y);

      gray_level = (double) (pixel.l * 255);

      pixel_.h = pixel.h;
      pixel.s = 0;
      pixel_.l = (double) C[gray_level]/N; //C[gray_level]*255/N;
      pixel_.a = pixel.a;
    }
  }

  return imgout;
}

//******************************************************************************************************************************************//
void updateNeighs(std::vector<Pt> &neighs, int i, int j) {
  neighs[0].setY(i-1); neighs[0].setX(j-1); 
  neighs[1].setY(i-1); neighs[1].setX(j); 
  neighs[2].setY(i-1); neighs[2].setX(j+1); 
  neighs[3].setY(i);   neighs[3].setX(j-1); 
  neighs[4].setY(i);   neighs[4].setX(j+1); 
  neighs[5].setY(i+1); neighs[5].setX(j-1); 
  neighs[6].setY(i+1); neighs[6].setX(j); 
  neighs[7].setY(i+1); neighs[7].setX(j+1); 
}

void updateCardNeighs(std::vector<Pt> &neighs, int i, int j) {
  neighs[0].setY(i-1);    neighs[0].setX(j);      //NORTH
  neighs[1].setY(i);      neighs[1].setX(j-1);    //WEST
  neighs[2].setY(i);      neighs[2].setX(j+1);    //EAST
  neighs[3].setY(i+1);    neighs[3].setX(j);      //SOUTH
}

bool isCard(Pt p, int i, int j) { //Is "p" a cardinal point for Pt(i,j) ?
  if(p.getY() == i && p.getX() == j-1) { //WEST
    return true;
  }
  if(p.getY() == i-1 && p.getX() == j) { //NORTH
    return true;
  }
  if(p.getY() == i && p.getX() == j+1) { //EAST
    return true;
  }
  if(p.getY() == i+1 && p.getX() == j) { //SOUTH
    return true;
  }
  return false;
}

bool belongsTo(std::vector<int> list, int n){
  for(unsigned int i=0; i<list.size(); i++){
    if(list[i] == n){
      return true;
    }
  }
  return false;
}

void LPE_updateNeigh(std::vector<std::vector<int>> &labels, std::vector<std::vector<int>> &stabilisation, int i, int j, int order){
  if(labels[i][j] == 0 || stabilisation[i][j] == 1){ //0 = WATERSHED
    return;
  }

  stabilisation[i][j] = 1; //FINAL S
  //Initialisation of neighbours list
  std::vector<Pt> neighs;
  Pt empty;
  for(int k = 0; k < 8; k++){
    neighs.push_back(empty);
  }    

  updateNeighs(neighs,i,j);

  for(unsigned k = 0; k < neighs.size(); k++){
    //Case : the neighbour does not exist or has already a label
    if( neighs[k].getY() < 0 || neighs[k].getX() < 0 || neighs[k].getY() >= (int)labels.size() || neighs[k].getX() >= (int)labels[0].size() || 
    stabilisation[neighs[k].getY()][neighs[k].getX()] != 0){ //0 = WATERSHED
      continue;
    }

    labels[neighs[k].getY()][neighs[k].getX()] = labels[i][j];
    stabilisation[neighs[k].getY()][neighs[k].getX()] = order;
  }
}

PNG LPE_homemade(PNG image, int tolerance) {
  int width = image.width(); 
  int height = image.height();

  uiuc::PNG tmp2;
  uiuc::PNG imgout(width,height);
  tmp2 = grayscale(image);
  tmp2 = egalisationHistogram(tmp2);

  //Initialisation of labels
  std::vector<std::vector<int>> labels(height, std::vector<int>(width, -1)); //INIT
  std::vector<std::vector<int>> stabilisation(height, std::vector<int>(width, 0)); //0 for "TO DO" & 1 for "DONE"

  //Initialisation of the first label on the first pixel
  int label = 1;
  labels[0][0] = label; 

  //Initialisation of neighbours list
  Pt tmp_neigh = Pt();
  int neigh_size = 8;
  std::vector<Pt> neighs(neigh_size,tmp_neigh);

  std::vector<int> allLabels;
  std::vector<std::vector<double>> labelsValues; //HSLA values associated to each label
  std::vector<double> hsla_values;
  std::vector<int> labelsValuesGray; //Gray values

  allLabels.push_back(0); //WATERSHED

  hsla_values.push_back(0); hsla_values.push_back(0); hsla_values.push_back(0);
  labelsValues.push_back(hsla_values); //BLACK COLOR
  hsla_values.clear();

  labelsValuesGray.push_back(0);

  allLabels.push_back(1); //FIRST LABEL

  HSLAPixel tmp_pixel = image.getPixel(0,0);
  hsla_values.push_back(tmp_pixel.h); 
  hsla_values.push_back(tmp_pixel.s); 
  hsla_values.push_back(tmp_pixel.l);
  labelsValues.push_back(hsla_values); //FIRST PIXEL VALUES
  hsla_values.clear();

  labelsValuesGray.push_back(tmp2.getPixel(0,0).h);

  //Labelise the entire image
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      updateNeighs(neighs,i,j);

      for(int k = 0; k < neigh_size; k++){
        //std::cout << "y " << neighs[k].getY() << " x " << neighs[k].getX() << std::endl;
        //Case : the neighbour does not exist or has already a label
        if( neighs[k].getY() < 0 || neighs[k].getX() < 0 || neighs[k].getY() >= height || neighs[k].getX() >= width || 
        labels[neighs[k].getY()][neighs[k].getX()] == 0 ){ //0 = WATERSHED
          continue;
        }
        //Same zone, same label
        if( labelsValuesGray[labels[i][j]] - tolerance <= tmp2.getPixel(neighs[k].getX(),neighs[k].getY()).h &&
          tmp2.getPixel(neighs[k].getX(),neighs[k].getY()).h <= labelsValuesGray[labels[i][j]] + tolerance &&
          labels[i][j] != 0 //WATERSHED
        ){
          labels[neighs[k].getY()][neighs[k].getX()] = labels[i][j];
          //std::cout << "labels[" << neighs[k].getY() << "][" << neighs[k].getX() << "] = " << labels[neighs[k].getY()][neighs[k].getX()] << std::endl;
          continue;
        }
        //If the center pixel is a WATERSHED then give a label for each neighbour
        if(labels[i][j] == 0){
          label++;
          labels[neighs[k].getY()][neighs[k].getX()] = label;

          //Determination of labels and its values
          allLabels.push_back(label); //label

          tmp_pixel = image.getPixel(neighs[k].getX(),neighs[k].getY());
          hsla_values.push_back(tmp_pixel.h); 
          hsla_values.push_back(tmp_pixel.s); 
          hsla_values.push_back(tmp_pixel.l);
          labelsValues.push_back(hsla_values); //DETERMINE A UNIQUE COLOR FOR THE LABEL
          hsla_values.clear();

          tmp_pixel = tmp2.getPixel(neighs[k].getX(),neighs[k].getY());
          labelsValuesGray.push_back(tmp_pixel.h);
          continue;
        }
        //Else it is a watershed
        labels[neighs[k].getY()][neighs[k].getX()] = 0;
      }
    }
  }
  //End of labelising

  /*for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      std::cout << "labels[" << i << "][" << j << "] = " << labels[i][j] << std::endl;
    }
  }*/
  
  //Stabilisation of the output
  int cpt = 1;

  //Stabilisation of the watersheds
  while(cpt != 0){
    cpt = 0;

    //Update and complete the watershed
    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){

        if(labels[i][j] == 0){ //Already a part of the watershed
          continue;
        }

        int nonExistingNeigh = 0;
        int nbNeighWatershed = 0;
        int nbCardNeighWatershed = 0;
        updateNeighs(neighs,i,j);

        for(int k = 0; k < neigh_size; k++){

          if( neighs[k].getY() >= 0 && neighs[k].getX() >= 0 && 
          neighs[k].getY() < height && neighs[k].getX() < width
          ){ 

            if(labels[neighs[k].getY()][neighs[k].getX()] == 0){ //0 = WATERSHED
              nbNeighWatershed++;

              if(isCard(neighs[k],i,j)){
                nbCardNeighWatershed++;
              }
            }
          }
          else{ //The neighbour does not exist
            nonExistingNeigh++;
          }

          if(nbNeighWatershed >= 6 || nbCardNeighWatershed >= 3 || (nbNeighWatershed >= 2 && nonExistingNeigh >= 3) ){ //The last : useful to handle minor errors on non-watershed pixels (which would be watersheded) 
            labels[i][j] = 0; //Become a part of the watershed
            cpt++;
            break;
          }
        }
      }
    }
  }

  //Initialisation of cardinal neighbours list
  std::vector<Pt> cardNeighs;
  for(int i = 0; i < 4; i++){
    cardNeighs.push_back(Pt());
  }

  //Stabilisation of the bassins
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      if(labels[i][j] == 0){ //WATERSHED
        stabilisation[i][j] = 1; //WATERSHED = STABILISED
      }
    }
  }

  cpt = 1;
  int cpt2 = 1;
  int order = 2;

  while(cpt != 0){ //cpt = nb of unstabilised + non-treated cases left
    cpt = 0;

    while(cpt2 != 0){ //cpt2 = nb of stabilised but next cases to be used left (considered as non-treated)
      cpt2 = 0;

      for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
          if(stabilisation[i][j] == order){ //order (case which will change)
            LPE_updateNeigh(labels,stabilisation,i,j,order+1); //Determination of next orders cases
            cpt2++;
            cpt++;
          }
        }
      }

      order++;
    }

    for(int i = 0; i < height; i++){
      if(cpt != 0){
        break;
      }

      for(int j = 0; j < width; j++){
        if(stabilisation[i][j] == 0){ //NOT STABILISED YET
          LPE_updateNeigh(labels,stabilisation,i,j,order); 
          cpt++;
          cpt2++; //To enter in the second while
          break;
        }
      }
    }
  }

  //Colorisation of the output image
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      HSLAPixel & pixel = imgout.getPixel(j,i);
      pixel.h = labelsValues[labels[i][j]][0];
      pixel.s = labelsValues[labels[i][j]][1];
      pixel.l = labelsValues[labels[i][j]][2];
      pixel.a = 1;
    }
  }

  return imgout;
}