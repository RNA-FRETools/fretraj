#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <stdlib.h>
#include <random>
#include <math.h> 

namespace py = pybind11;

std::random_device rd;
std::mt19937 gen(rd());


int findExcitationIndex(int traj_length, int skipframesatstart, int skipframesatend){
    std::uniform_int_distribution<> excitationIndex(skipframesatstart, traj_length-skipframesatend);
    return excitationIndex(gen);
}

int integerChoice(std::vector<double> probabilities){
    std::discrete_distribution<> sample(probabilities.begin(), probabilities.end());
    return sample(gen);
}

int polarization(std::vector<double> excitation_dipole, std::vector<double> emission_dipole){
    double product = std::inner_product(excitation_dipole.begin(), excitation_dipole.end(), emission_dipole.begin(), 0.0);
    double product_square = product*product;
    double prob_parallel = 2*product_square/(product_square+1);
    std::vector<double> probabilities {prob_parallel, 1-prob_parallel};
    std::discrete_distribution<> pol(probabilities.begin(), probabilities.end());
    return pol(gen);
}

std::tuple<int,int> findRelaxationIndex_DD_DA(double pD_tot, py::array_t<double> pD_totfret, int excitation_ndx, int traj_length){
    auto pD_totfret_uc = pD_totfret.unchecked<1>();
    for(int ndx=excitation_ndx; ndx<traj_length; ++ndx) {
        double p = std::generate_canonical<double, 10>(gen);
        if (p < pD_tot){
            return std::make_tuple(1, ndx);
        }
        if (p < pD_totfret_uc(ndx)){
            return std::make_tuple(2, ndx);
        }
    }
    return std::make_tuple(0, traj_length);
}

std::tuple<int,int> findRelaxationIndex_AA(double pA_tot, int excitation_ndx, int traj_length){
    for(int ndx=excitation_ndx; ndx<traj_length; ++ndx) {
        double p = std::generate_canonical<double, 10>(gen);
        if (p < pA_tot){
            return std::make_tuple(2, ndx);
        }
    }
    return std::make_tuple(0, traj_length);
}

int checkRelaxationIndex(int event, int excitation_ndx, int relaxation_ndx, py::array_t<double> R, double quenching_radius, double QD, double QA){
    double p = std::generate_canonical<double, 10>(gen);
    auto r = R.unchecked<1>();
    for (int i = excitation_ndx; i <= relaxation_ndx; i++){
        if (r(i) < quenching_radius){
            return 0;
        }
    }
    if (event == 1) {
        if (p > QD){
            return -1;
        }
    } else {
        if (p > QA){
            return -2;
        }
    }
    return event;
}


PYBIND11_MODULE(relaxation, m) {
    m.doc() = "C-extension of the fluordynamics package for fast photon sampling";
    m.def("findRelaxationIndex_DD_DA", &findRelaxationIndex_DD_DA, "Return a relaxation event upon donor excitation from a donor dye (1, photon/IC), an acceptor dye (1, photon/IC) or no relaxation (0, remain in excited state)", 
        py::arg("pD_tot"), py::arg("pD_totfret"), py::arg("excitation_ndx"), py::arg("traj_length"));
    m.def("findRelaxationIndex_AA", &findRelaxationIndex_AA, "Return a relaxation event upon acceptor excitation from an acceptor dye (2, photon/IC) or no relaxation (0, remain in excited state)", 
        py::arg("pD_tot"), py::arg("excitation_ndx"), py::arg("traj_length"));
    m.def("checkRelaxationIndex", &checkRelaxationIndex, "Check if the relaxation event is due to quenching (0), emission of a donor photon (1), an acceptor photon (2) or "
                                                        "internal conversion from the donor (-1) of acceptor (-2) dye. For this purpose, the relaxation event is first evaluated "
                                                        "with respect to the inter-dye distance. If at any point in the interval between excitation_ndx and relaxation_ndx the distance "
                                                        "is smaller than the quenching radius the event is defined as an internal conversion (i.e. no photon emission). " 
                                                        "Next, the probability of a donor or acceptor photon is checked against the quantum yield of the dyes.)", 
                                                        py::arg("event"), py::arg("excitation_ndx"), py::arg("relaxation_ndx"), py::arg("R"), py::arg("quenching_radius"), py::arg("QD"), py::arg("QA"));
    m.def("findExcitationIndex", &findExcitationIndex, "Return an excitation event within the trajectory interval defined by \"skipframesatstart\" and \"skipframesatend\"", 
        py::arg("traj_length"), py::arg("skipframesatstart"), py::arg("skipframesatend"));
    m.def("integerChoice", &integerChoice, "Return an integer in the interval [0,n) with the specified discrete probabilities [p_1,p_2,...,p_n]", py::arg("probabilities"));
    m.def("polarization", &polarization, "Return the polarization of the photon where 0 is a parallel p-photon and 1 is a perpendicular s-photon"
                                         "The intensity of the parallel polarized light is proportional to cos^2(x), where x is the angle between between the dipoles at time points t_ex and t_em."
                                         "Similarly, the intensity of perpedicular polarized light is proportional to 0.5*sin^2(x)"
                                         "The probability of a parallel photon is thus calculated as cos^2(x)/(cos^2(x)+0.5*sin^2(x)) = 2*cos^2(x)/(cos^2(x)+1)"),
                                         py::arg("excitation_dipole"), py::arg("emission_dipole");
}