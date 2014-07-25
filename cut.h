using namespace std;
const int BKLNUMBER =30;
const int MODULENUMER =5;
const int TacdiffEastUp[BKLNUMBER][MODULENUMER] = {{2000}};
const int TacdiffEastbottom[BKLNUMBER][MODULENUMER] = {{2000}};
const int TacdiffWestUp[BKLNUMBER][MODULENUMER] = {{2000}};
const int TacdiffWestbottom[BKLNUMBER][MODULENUMER] = {{2000}};
const float PT = 1.0;
const float ETA[2] = {-0.8,0.8};
const int NFTPTS = 15;
const int NDEDXPTS = 10;
const float NSIGMAPI[2] ={-1.0,3.0};
const float mvDrift =60.;
const int nNeighborsFlag = 0;//0 for kFALSE; 1 for kTRUE
