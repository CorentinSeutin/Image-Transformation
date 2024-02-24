#include "uiuc/PNG.h"
using namespace uiuc;

class Pt{
    public:
        int x, y;

        Pt();
        Pt(int j, int i);
        
        int getX();
        int getY();
        void setX(int j);
        void setY(int i);
};