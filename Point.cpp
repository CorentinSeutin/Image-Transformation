#include "Point.hpp"

Point::Point() {
    x = 0; y = 0;
}

Point::Point(unsigned j, unsigned i) {
    x = j; y = i;
}

unsigned Point::getX() { return x; }
unsigned Point::getY() { return y; }
void Point::setX(unsigned j){ x = j; }
void Point::setY(unsigned i){ y = i; }