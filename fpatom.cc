#include <madness/world/world.h>
#include <madness/mra/mra.h>
#include <apps/chem/xcfunctional.h>

using namespace madness;

class NuclearPotentialFunctor: public FunctionFunctorInterface<double, 3> 
{
    private:

        int zn;

    public:

        NuclearPotentialFunctor(int zn) : zn(zn)
        {
        }

        double operator()(const coord_3d& x) const 
        {
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            double r0 = 1e-8;
            if (r > r0) 
            {
                return -zn / r;
            }
            else
            {
                double a = 0.5 * zn / std::pow(r0, 3);
                double c = -1.5 * zn / r0; 
                return a * r * r + c;
            }
        }
};

class AtomicOrbitalFunctor: public FunctionFunctorInterface<double, 3> 
{
    private:
        
        int zn;

        int n;

        int l;

        int m;

    double LaguerreL(double x) const
    {
        return 1.0;
    }

    public:

        AtomicOrbitalFunctor(int zn, int n, int l, int m) : zn(zn), n(n), l(l), m(m)
        {
        }

        double operator()(const coord_3d& x) const 
        {
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

            return std::pow(2 * zn * r / n, l) * std::exp(-zn * r / n) * LaguerreL(2 * zn * r / n);
        }
};

static double zero_3d(const coord_3d& r)
{
    return 0.0;
}

void test_xc(World& world)
{
    bool spin_polarized = true;

    XCfunctional xcfunc;

    xcfunc.initialize("LDA", spin_polarized, world);

    int npt = 5;
    double rho_up[] = {0, 0.1, 0.2, 0.5, 2.234};
    double rho_dn[] = {0, 0.0, 0.1, 1.5, 2.234};

    Tensor<double> rhoa_t(npt);
    Tensor<double> rhob_t(npt);

    for (int i = 0; i < npt; i++)
    {
        rhoa_t(i) = rho_up[i];
        rhob_t(i) = rho_dn[i];
    }

    std::vector< Tensor<double> > xc_args;
    xc_args.push_back(rhoa_t);
    xc_args.push_back(rhob_t);

    Tensor<double> vxc_up;
    Tensor<double> vxc_dn;
    Tensor<double> exc;

    exc = xcfunc.exc(xc_args, -1);
    vxc_up = xcfunc.vxc(xc_args, 0, 0);
    vxc_dn = xcfunc.vxc(xc_args, 1, 0);

    for (int i = 0; i < npt; i++)
    {
        printf("%12.6f %12.6f %12.6f %12.6f %12.6f\n", rhoa_t(i), rhob_t(i), vxc_up(i), vxc_dn(i), exc(i));
    }
}

//real_function_3d make_dft_potential(const XCfunctional& xc, const vector_real_function_3d& vf, int ispin, int what)
//{
//  return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
//}

double make_dft_energy(const XCfunctional& xc, const vector_real_function_3d& rho)
{
    /* spin is not used in xc_functional, so we can pass any value */
    real_function_3d exc_rho = multiop_values<double, xc_functional, 3>(xc_functional(xc, -1), rho);
    return exc_rho.trace();
}

int main(int argc, char** argv) 
{
    auto& world = initialize(argc, argv);
    /////World world(SafeMPI::COMM_WORLD);
    
    startup(world, argc, argv);

    // MRA parameters
    double thresh = 1e-4;
    int kwavelet = 6;
    int truncate_mode = 0; 
    double L = 50.0;

    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_k(kwavelet);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_truncate_mode(truncate_mode);

    /* init spin-polarized XC potential */
    XCfunctional xcfunc;
    xcfunc.initialize("LDA", true, world);

    real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(new NuclearPotentialFunctor(1))).thresh(1e-10).truncate_on_project();
    //vnuc.reconstruct();

    real_function_3d psi_i = real_factory_3d(world).functor(real_functor_3d(new AtomicOrbitalFunctor(1, 1, 0, 0))).thresh(1e-10).truncate_on_project();
    psi_i.scale(1.0/psi_i.norm2());


    vector_real_function_3d rho(2);
    rho[0] = real_factory_3d(world).f(zero_3d);
    rho[1] = real_factory_3d(world).f(zero_3d);
    
    rho[0] += (psi_i * psi_i);
    rho[1] += (psi_i * psi_i); // If I comment this, the code crashes. Why?

    printf("total electron charge: %18.12f\n", rho[0].trace() + rho[1].trace());

    double Exc = make_dft_energy(xcfunc, rho);

    printf("Exc: %18.12f\n", Exc);


    //double Ekin = 0;
    //std::vector< std::shared_ptr<real_derivative_3d> > gradop = gradient_operator<double, 3>(world);
    //for (int axis = 0; axis < 3; axis++)
    //{
    //    real_function_3d dpsi = apply(*(gradop[axis]), psi_i, true);
    //    Ekin += inner(dpsi, dpsi);
    //}
    //Ekin *= 0.5;

    //double Epot = inner(rho, vnuc);

    //printf("Total energy: %18.12f\n", Ekin + Epot);
    //
    //return 0;
}
