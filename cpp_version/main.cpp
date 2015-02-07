#include "EVAP.h"

using namespace std;

int main() {
    double lat = 51.4;
    int n = 172;
    double elv = 74.0; 
    int y = 2001; 
    float sf = 0.43; 
    float tc = 17.3; 
    float sw = 0.5;
    EVAP my_evap(lat, n, elv, y, sf, tc, sw);
    my_evap.display();
}
