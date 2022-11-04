#pragma once

namespace myLib
{
    const double calcHopfQ(const double &tau);
    const double calcTemperature(const double &tau,
                                 const double &Teff);
    const double planck(const double &lam,
                        const double &T);
}