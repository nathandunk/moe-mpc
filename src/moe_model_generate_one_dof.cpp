#include <Mahi/Casadi/ModelGenerator.hpp>
#include <Mahi/Casadi/M.hpp>
#include <Mahi/Casadi/G.hpp>
#include <Mahi/Casadi/V.hpp>
#include <Mahi/Util.hpp>

using namespace casadi;
using mahi::util::PI;
using mahi::util::DEG2RAD;

int main(int argc, char* argv[])
{
    mahi::util::Options options("options.exe", "Simple Program Demonstrating Options");

    options.add_options()
        ("l,linear", "Generates linearized model.")
        ("j,joint_num", "joint number to generate model for (0 indexed)", mahi::util::value<int>())
        ("h,help", "Prints help message.");
    
    auto result = options.parse(argc, argv);

    if (result.count("help")){
        std::cout << options.help() << std::endl;
        return 0;
    }

    bool linear = result.count("linear") > 0;
    int joint_num = result.count("joint_num") ? result["joint_num"].as<int>() : 0;

    SX x, x_dot, u;
    std::string model_name;

    model_name = "moe_j" + std::to_string(joint_num);
    SX q0 = SX::sym("q0");
    SX q1 = SX::sym("q1");
    SX q2 = SX::sym("q2");
    SX q3 = SX::sym("q3");
    SX qs = SX::vertcat({q0,q1,q2,q3});
    SX q  = qs(joint_num);


    SX q0_dot = SX::sym("q0_dot");
    SX q1_dot = SX::sym("q1_dot");
    SX q2_dot = SX::sym("q2_dot");
    SX q3_dot = SX::sym("q3_dot");
    SX q_dots = SX::vertcat({q0_dot,q1_dot,q2_dot,q3_dot});
    SX q_dot  = q_dots(joint_num);
    
    SX T0 = SX::sym("T0");
    SX T1 = SX::sym("T1");
    SX T2 = SX::sym("T2");
    SX T3 = SX::sym("T3");
    SX Ts  = SX::vertcat({T0, T1, T2, T3});
    SX T = Ts(joint_num);

    // inertia from each joint to the rest of the robot when all DOFs are at 0 (kg*m^2)
    std::vector<double> I_rests = {0.18509398, // joint 0
                                   1.0, 
                                   1.0, 
                                   1.0};

    SX I = I_rests[joint_num];

    // state vector
    x = SX::vertcat({q, q_dot});
    u = SX::vertcat({T});

    std::vector<double> qes = {+33.7023*DEG2RAD, 0.0, 0.0, 0.0}; // rad 
    std::vector<double> rs  = {0.0286, 0.0, 0.0, 0.0}; // m
    std::vector<double> ms  = {9.55877228, 0.0, 0.0, 0.0}; // kg

    SX qe = qes[joint_num];
    SX r = rs[joint_num];
    SX m = ms[joint_num];

    std::vector<double>  B_coef = {0.0840, 0.0252, 0.0019, 0.0029};
    std::vector<double> Fk_coef = {0.2186, 0.1891, 0.0541, 0.1339};

    SX B = B_coef[joint_num]*q_dot(joint_num)*1.0;
    SX Fk = Fk_coef[joint_num]*tanh(q_dot(joint_num)*10.0);

    // Transmission Rations [rad/m]
    std::vector<double> eta = {0.4200/4.50, 0.4706/8.75, 0.4735/9.00, 0.2210/6.00};
    // Motor rotor inertias [Kg*m^2]
    std::vector<double> Jm_ = {1340e-7, 137e-7, 137e-7, 34.7e-7};

    SX Jm = Jm_[joint_num]/eta[joint_num]/eta[joint_num];

    SX G = r*m*g*cos(15.0*DEG2RAD)*sin(q+qe);

    SX q_d_dot = (T - (B+Fk+G))/(Jm+I);
    x_dot = vertcat(q_dot,q_d_dot);

    if (linear) model_name = "linear_" + model_name;
    
    // Bounds on state
    std::vector<double> x_min(x.size1(),-inf);
    std::vector<double> x_max(x.size1(), inf);

    // Bounds for control
    std::vector<double> u_min(u.size1(),-inf);
    std::vector<double> u_max(u.size1(), inf);

    // settings for multiple shooting constructions
    mahi::util::Time time_step  = mahi::util::milliseconds(linear ? 2 : 10);
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
