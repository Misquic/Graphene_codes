#include "funkcje.h"
#include "Constants.h" // freqquently and widely used variables used as parameters



Bilayer::Bilayer(double nit, double nib,
                 double dt, double dg, double db,
                 double et, double eg, double eb):
    m_nit{nit*Const::inv_cmsq2au}, m_nib{nib*Const::inv_cmsq2au},
    m_dt{std::abs(dt)*Const::nm2au}, m_dg{std::abs(dg)*Const::nm2au}, m_db{std::abs(db)*Const::nm2au},
    m_et{std::abs(et)}, m_eg{std::abs(eg)}, m_eb{std::abs(eb)},
    m_Ct{Const::eps_0*m_et/m_dt}, m_Cg{Const::eps_0*m_eg/m_dg}, m_Cb{Const::eps_0*m_eb/m_db} 
{
};

inline void count_n_V(double V, double Vg_o, double C, double m_Cg, double abs_ni, double sign_ni,
                      double &n, double &Vg)
{
    double nQ      = M_PI_2 * pow2(Const::hv * (C + m_Cg) );
    double sqrt_2_nQ_abs_ni = std::sqrt(2*nQ*abs_ni);
    double nC      = sign_ni * abs_ni + (C*V + m_Cg*Vg_o) + sign_ni * sqrt_2_nQ_abs_ni;
    double delta_n = signum(nC) * nQ * ( 1 - std::sqrt(1 + 2*std::abs(nC)/nQ));
    
    n = nC + delta_n;
    Vg = - (delta_n + sign_ni*sqrt_2_nQ_abs_ni) / (C + m_Cg);
}

std::pair<double, double> Bilayer::count_densities(double Vt, double Vb)
{
    constexpr double tol = 1e-4;
    Vt *= Const::eV2au;
    Vb *= Const::eV2au;
    
    double Vgt     = 0, Vgb     = 0;
    double Vgt_old = 0, Vgb_old = 0;
    double nt      = 0, nb      = 0;
    const double abs_nit = std::abs(m_nit);
    const double abs_nib = std::abs(m_nib);
    const double sign_nit = m_nit >= 0. ? 1. : -1.;
    const double sign_nib = m_nib >= 0. ? 1. : -1.;

    unsigned it = 0;
    bool converged = false;
    do
    {
        it++;
        // Vgt = Vgt_old;
        // Vgb = Vgb_old;
        // t
        count_n_V(Vt, Vgb_old, m_Ct, m_Cg, abs_nit, sign_nit, nt, Vgt);
        // b
        count_n_V(Vb, Vgt_old, m_Cb, m_Cg, abs_nib, sign_nib, nb, Vgb);
        
        
        if(it % 15 == 0)
        {
            if ((Vgt - Vgt_old)/Vgt_old < tol && (Vgb - Vgb_old)/Vgb_old < tol)
            {
                converged = true;
            }
            // printf("\r         \r%i", it);
        }
        Vgt_old = Vgt;
        Vgb_old = Vgb;
    } while ( !converged );

    // printf("iterations: %i\n", it);

    return {nt, nb};
};


//
