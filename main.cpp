/**
 * @file main.cpp
 * A simple C++ program that manipulates an image.
 *
 * @author University of Illinois CS 225 Course Staff
 * @author Updated by University of Illinois CS 400 Course Staff
**/

#include "ImageTransform.h"
#include "uiuc/PNG.h"

int main() {
  uiuc::PNG png, result; //png2;
  //uiuc::PNG image(2480,3508);

  //png.readFromFile("alma.png");
  png.readFromFile("radio.png");
  //png2.readFromFile("accolade.png");
  
  /*
  result = grayscale(png);
  result.writeToFile("out-grayscale.png");
  
  result = createSpotlight(png, 450, 150);
  result.writeToFile("out-spotlight.png");

  result = illinify(png);
  result.writeToFile("out-illinify.png");*/

  //png2.readFromFile("overlay.png");
  //result = watermark(png, png2);
  //result.writeToFile("out-watermark.png");

  //result = HFilter(png);
  //result.writeToFile("partition.png");

  //image = Partition(image);
  //image.writeToFile("partition.png");

  //result = Divide(png2);
  //result = DivideByThree(result);
  //result = RemoveWhite(result);
  //result.writeToFile("CurlyBracket.png");

  //result = logoEUREKA6(png);
  //result.writeToFile("eureka_res.png");

  result = LPE_homemade(png,20);
  //result = egalisationHistogram(png);
  result.writeToFile("AAA.png");

  /*try{
    result = createRandomGraph(1000000,1333,1333);
    result.writeToFile("out-createRandomGraph3x3.png");
  }
  catch(std::invalid_argument& e){
    std::cerr << e.what() << std::endl;
    return -1;
  }*/

  return 0;
}
