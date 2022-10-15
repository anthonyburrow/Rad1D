#include <iostream>

#include "gaussianQuadrature.hpp"
#include "RadModel.hpp"

using namespace std;

namespace myLib
{
    void getWeights(myLib::RadModel &radModel)
    {
        const int &nQuad = radModel.params.nQuad;

        switch(radModel.params.nQuad)
        {
            case 2:
            {
                const double inv_sqrt_3 = 1.0 / sqrt(3.0);
                radModel.quadMu = { inv_sqrt_3, -inv_sqrt_3 };
                radModel.quadW = { 1.0, 1.0 };
                break;
            }
            case 4:
            {
                radModel.quadMu = { 0.8611363115940526, 0.3399810435848563, -0.3399810435848563, -0.8611363115940526 };
                radModel.quadW = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
                break;
            }
            case 6:
            {
                radModel.quadMu = { 0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.2386191860831969,
                                    -0.6612093864662645, -0.9324695142031521 };
                radModel.quadW = { 0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910,
                                   0.3607615730481386, 0.1713244923791704 };
                break;
            }
            case 8:
            {
                radModel.quadMu = { 0.9602898564975363, 0.7966664774136267, 0.5255324099163290, 0.1834346424956498,
                                    -0.1834346424956498, -0.5255324099163290, -0.7966664774136267, -0.9602898564975363 };
                radModel.quadW = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                   0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
                break;
            }
            default:
                cout << "nQuad = " << nQuad << " is not valid" << endl;
        }
    }
}