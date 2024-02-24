#include "uiuc/PNG.h"
using namespace uiuc;

class Point{
    public:
        unsigned x, y;

        Point();
        Point(unsigned j, unsigned i);
        
        unsigned getX();
        unsigned getY();
        void setX(unsigned j);
        void setY(unsigned i);
};