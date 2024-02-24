#include "Pt.hpp"

Pt::Pt() {
    x = 0; y = 0;
}

Pt::Pt(int j, int i) {
    x = j; y = i;
}

int Pt::getX() { return x; }
int Pt::getY() { return y; }
void Pt::setX(int j){ x = j; }
void Pt::setY(int i){ y = i; }