#include <Mahi/Casadi/ModelGenerator.hpp>
#include <Mahi/Casadi/M.hpp>
#include <Mahi/Casadi/G.hpp>
#include <Mahi/Casadi/V.hpp>
#include <Mahi/Util.hpp>

using namespace casadi;
using mahi::util::PI;

int main(int argc, char* argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        ("l,linear", "Generates linearized model.")
        ("h,help", "Prints help message.");
    
    auto result = options.parse(argc, argv);

    if (result.count("help")){
        std::cout << options.help() << std::endl;
        return 0;
    }

    bool linear = result.count("linear") > 0;

    SX x, x_dot, u;
    std::string model_name;

    model_name = "moe";
    SX q0 = SX::sym("q0");
    SX q1 = SX::sym("q1");
    SX q2 = SX::sym("q2");
    SX q3 = SX::sym("q3");
    SX q0_dot = SX::sym("q0_dot");
    SX q1_dot = SX::sym("q1_dot");
    SX q2_dot = SX::sym("q2_dot");
    SX q3_dot = SX::sym("q3_dot");
    
    SX T0 = SX::sym("T0");
    SX T1 = SX::sym("T1");
    SX T2 = SX::sym("T2");
    SX T3 = SX::sym("T3");

    // state vector
    x = SX::vertcat({q0,q1,q2,q3,q0_dot,q1_dot,q2_dot,q3_dot});
    u = SX::vertcat({T0,T1,T2,T3});

    SX q_dot = SX::vertcat({q0_dot,q1_dot,q2_dot,q3_dot});

    std::vector<double>  B_coef = {0.1215, 0.0252, 0.0019, 0.0029};
    std::vector<double> Fk_coef = {   0.5, 0.1891, 0.0541, 0.1339};

    SX B = SX::vertcat({B_coef[0]*q0_dot*1.0, B_coef[1]*q1_dot*1.0, B_coef[2]*q2_dot*1.0, B_coef[3]*q3_dot*1.0});
    SX Fk = SX::vertcat({Fk_coef[0]*tanh(q0_dot*10.0), Fk_coef[1]*tanh(q1_dot*10.0), Fk_coef[2]*tanh(q2_dot*10.0), Fk_coef[3]*tanh(q3_dot*10.0)});

    // Fk[0] = Fk_coef[0]*std::tanh(q0d*10);
    // Fk[1] = Fk_coef[1]*std::tanh(q1d*10);
    // Fk[2] = Fk_coef[2]*std::tanh(q2d*10);
    // Fk[3] = Fk_coef[3]*std::tanh(q3d*10);

    auto V = get_V(x);
    auto G = get_G(x);
    auto B_eom = u - mtimes(V,q_dot) - G - B - Fk;
    auto A_eom = get_M(x);
    SX q_d_dot = solve(A_eom,B_eom);
    x_dot = vertcat(q_dot,q_d_dot);

    if (linear) model_name = "linear_" + model_name;
    
    // Bounds on state
    std::vector<double> x_min(x.size1(),-inf);
    std::vector<double> x_max(x.size1(), inf);

    // Bounds for control
    std::vector<double> u_min(u.size1(),-inf);
    std::vector<double> u_max(u.size1(), inf);

    // settings for multiple shooting constructions
    mahi::util::Time time_step  = mahi::util::milliseconds(2);
    int num_shooting_nodes = 25;

    ModelParameters model_parameters(model_name, // name
                                     x.size1(),                // num_x
                                     u.size1(),                // num_u
                                     time_step,                // step_size
                                     num_shooting_nodes,       // num_shooting_nodes
                                     linear);                  // is_linear;                  

    // 
    ModelGenerator my_generator(model_parameters, x, x_dot, u);

    my_generator.create_model();
    my_generator.generate_c_code();
    my_generator.compile_model();

    return 0;
}
