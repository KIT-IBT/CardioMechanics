/* -------------------------------------------------------
 
 CBTensionModelLand17.cpp
 
 Ver. 1.0.0
 
 Created:       Tobias Gerach  (24.07.2019)
 Last modified: Tobias Gerach  (08.01.2023)
 
 Institute of Biomedical Engineering
 Karlsruhe Institute of Technology (KIT)
 
 http://www.ibt.kit.edu
 
 Copyright 2000-2009 - All rights reserved.
 
 ------------------------------------------------------ */

#include "CBTensionModelLand17.h"
#include <algorithm>
#include <cassert>

TFloat CBTensionModelLand17::dXSdt(TFloat xb_ws, TFloat xb_su, TFloat xb_su_gamma) {
    return xb_ws - xb_su - xb_su_gamma;
}

TFloat CBTensionModelLand17::dXWdt(TFloat xb_uw, TFloat xb_wu, TFloat xb_ws, TFloat xb_wu_gamma) {
    return xb_uw - xb_wu - xb_ws - xb_wu_gamma;
}

TFloat CBTensionModelLand17::dTRPNdt(TFloat TRPNk, TFloat Cai, TFloat Ca50, TFloat TRPNn, TFloat TRPN) {
    return TRPNk * (pow((Cai/Ca50), TRPNn) * (1.0 - TRPN) - TRPN);
}

TFloat CBTensionModelLand17::dTmBlockeddt(TFloat ktm_block, TFloat TRPN_NP, TFloat XU, TFloat ku, TFloat TRPN,
                                          TFloat nTM, TFloat TmBlocked) {
    return ktm_block * TRPN_NP * XU - ku * pow(TRPN, nTM/2.0) * TmBlocked;
}

TFloat CBTensionModelLand17::dZetaSdt(TFloat A, TFloat dlambdadt, TFloat cs, TFloat ZetaS) {
    return A * dlambdadt - cs * ZetaS;
}

TFloat CBTensionModelLand17::dZetaWdt(TFloat A, TFloat dlambdadt, TFloat cw, TFloat ZetaW) {
    return A * dlambdadt - cw * ZetaW;
}

TFloat CBTensionModelLand17::dC_ddt(TFloat k, TFloat eta, TFloat C_s) {
    return k * C_s / eta;
}

TFloat CBTensionModelLand17::CalcActiveStiffness() {
    TFloat lambda_m = std::min(1.2, S_.lambda);
    TFloat h = std::max(0.0, Overlap(lambda_m));
    
    return h * Tref_/rs_ * (A_ * S_.XS + A_ * S_.XW);
}

TFloat CBTensionModelLand17::CaBrixius(TFloat t) {
    // Calcium transient that was used in the Land18 publication.
    // Values were extracted from the paper and interpolated to ms time-steps.
    t = (t - GetActivationTime());
    while (t > cycleLength_) {
        t = t - cycleLength_;
    }
    int index = t/1e-3;
    if (index < 0) {
        index = 0;
    }
    
    // if cycleLength > 1s, use last index
    index = std::min(1000, index);
    double humanCai[] = {
        0.100000, 0.136592, 0.169485, 0.228262, 0.256097, 0.268791, 0.318309, 0.386736, 0.417231, 0.445595, 0.461671,
        0.475009, 0.482842, 0.490374, 0.494654, 0.498933, 0.503213, 0.507444, 0.511293, 0.515142, 0.518992, 0.523474,
        0.528889, 0.534304, 0.539701, 0.545097, 0.549805, 0.551967, 0.554129, 0.556291, 0.558454, 0.560616, 0.562857,
        0.566035, 0.569213, 0.572391, 0.575569, 0.578747, 0.580728, 0.582555, 0.584382, 0.586208, 0.588035, 0.589862,
        0.591689, 0.592770, 0.593671, 0.594573, 0.595475, 0.596376, 0.597278, 0.598179, 0.598724, 0.598225, 0.597726,
        0.597226, 0.596727, 0.596228, 0.595728, 0.595229, 0.594454, 0.593561, 0.592669, 0.591777, 0.590884, 0.589992,
        0.589100, 0.588213, 0.587358, 0.586502, 0.585647, 0.584792, 0.583936, 0.583081, 0.582226, 0.581210, 0.580086,
        0.578963, 0.577839, 0.576716, 0.575592, 0.574469, 0.573355, 0.572440, 0.571525, 0.570609, 0.569694, 0.568778,
        0.567863, 0.566947, 0.565860, 0.564598, 0.563335, 0.562073, 0.560811, 0.559549, 0.558287, 0.557025, 0.555671,
        0.554311, 0.552952, 0.551593, 0.550234, 0.548874, 0.547515, 0.546159, 0.544806, 0.543454, 0.542102, 0.540749,
        0.539397, 0.538045, 0.536692, 0.535270, 0.533834, 0.532399, 0.530963, 0.529528, 0.528092, 0.526657, 0.525208,
        0.523730, 0.522251, 0.520773, 0.519294, 0.517816, 0.516337, 0.514859, 0.513269, 0.511641, 0.510013, 0.508384,
        0.506756, 0.505128, 0.503499, 0.501877, 0.500282, 0.498687, 0.497092, 0.495497, 0.493902, 0.492307, 0.490712,
        0.488993, 0.487204, 0.485415, 0.483626, 0.481836, 0.480047, 0.478258, 0.476469, 0.474679, 0.472890, 0.471101,
        0.469312, 0.467523, 0.465733, 0.463944, 0.462106, 0.460225, 0.458344, 0.456463, 0.454582, 0.452701, 0.450820,
        0.448941, 0.447064, 0.445187, 0.443310, 0.441433, 0.439557, 0.437680, 0.435842, 0.434069, 0.432297, 0.430524,
        0.428752, 0.426979, 0.425206, 0.423434, 0.421736, 0.420055, 0.418374, 0.416693, 0.415012, 0.413331, 0.411650,
        0.410030, 0.408574, 0.407117, 0.405661, 0.404205, 0.402748, 0.401292, 0.399836, 0.398469, 0.397137, 0.395806,
        0.394474, 0.393143, 0.391811, 0.390480, 0.389141, 0.387768, 0.386395, 0.385022, 0.383649, 0.382276, 0.380903,
        0.379529, 0.378139, 0.376738, 0.375337, 0.373937, 0.372536, 0.371135, 0.369734, 0.368344, 0.367119, 0.365894,
        0.364669, 0.363444, 0.362219, 0.360993, 0.359768, 0.358536, 0.357296, 0.356057, 0.354817, 0.353578, 0.352338,
        0.351099, 0.349859, 0.348527, 0.347191, 0.345856, 0.344520, 0.343184, 0.341849, 0.340513, 0.339321, 0.338336,
        0.337352, 0.336367, 0.335382, 0.334397, 0.333413, 0.332428, 0.331487, 0.330553, 0.329619, 0.328685, 0.327751,
        0.326817, 0.325883, 0.324863, 0.323649, 0.322435, 0.321221, 0.320007, 0.318793, 0.317579, 0.316365, 0.315293,
        0.314267, 0.313241, 0.312216, 0.311190, 0.310164, 0.309138, 0.308087, 0.306935, 0.305783, 0.304630, 0.303478,
        0.302326, 0.301174, 0.300022, 0.298883, 0.297745, 0.296607, 0.295469, 0.294332, 0.293285, 0.292390, 0.291494,
        0.290599, 0.289703, 0.288808, 0.287698, 0.286460, 0.285221, 0.283983, 0.282745, 0.281507, 0.280268, 0.279252,
        0.278322, 0.277393, 0.276464, 0.275535, 0.274605, 0.273676, 0.272758, 0.271898, 0.271038, 0.270178, 0.269318,
        0.268459, 0.267599, 0.266739, 0.265828, 0.264884, 0.263941, 0.262998, 0.262055, 0.261112, 0.260169, 0.259235,
        0.258444, 0.257654, 0.256863, 0.256073, 0.255282, 0.254492, 0.253701, 0.252981, 0.252329, 0.251677, 0.251026,
        0.250374, 0.249722, 0.249070, 0.248418, 0.247660, 0.246897, 0.246134, 0.245371, 0.244608, 0.243845, 0.243083,
        0.242337, 0.241615, 0.240894, 0.240173, 0.239452, 0.238730, 0.238009, 0.237288, 0.236682, 0.236094, 0.235507,
        0.234920, 0.234333, 0.233746, 0.233159, 0.232557, 0.231924, 0.231291, 0.230657, 0.230024, 0.229390, 0.228757,
        0.228124, 0.227494, 0.226865, 0.226236, 0.225607, 0.224979, 0.224350, 0.223721, 0.223091, 0.222458, 0.221825,
        0.221191, 0.220558, 0.219924, 0.219291, 0.218658, 0.217982, 0.217284, 0.216586, 0.215888, 0.215190, 0.214491,
        0.213793, 0.213104, 0.212499, 0.211893, 0.211287, 0.210682, 0.210076, 0.209470, 0.208865, 0.208289, 0.207739,
        0.207189, 0.206639, 0.206089, 0.205539, 0.204988, 0.204438, 0.203874, 0.203310, 0.202746, 0.202182, 0.201618,
        0.201054, 0.200490, 0.199938, 0.199402, 0.198866, 0.198329, 0.197793, 0.197257, 0.196720, 0.196184, 0.195635,
        0.195085, 0.194535, 0.193985, 0.193435, 0.192884, 0.192334, 0.191819, 0.191370, 0.190922, 0.190474, 0.190025,
        0.189577, 0.189128, 0.188680, 0.188213, 0.187741, 0.187270, 0.186798, 0.186327, 0.185855, 0.185383, 0.184915,
        0.184457, 0.184000, 0.183542, 0.183084, 0.182627, 0.182169, 0.181711, 0.181282, 0.180866, 0.180450, 0.180034,
        0.179618, 0.179202, 0.178785, 0.178366, 0.177922, 0.177478, 0.177034, 0.176590, 0.176146, 0.175703, 0.175259,
        0.174831, 0.174415, 0.173999, 0.173583, 0.173167, 0.172751, 0.172334, 0.171919, 0.171512, 0.171105, 0.170698,
        0.170291, 0.169884, 0.169478, 0.169071, 0.168686, 0.168325, 0.167965, 0.167604, 0.167244, 0.166883, 0.166522,
        0.166162, 0.165767, 0.165369, 0.164972, 0.164574, 0.164176, 0.163779, 0.163381, 0.163024, 0.162732, 0.162441,
        0.162150, 0.161859, 0.161567, 0.161276, 0.160985, 0.160636, 0.160276, 0.159915, 0.159554, 0.159194, 0.158833,
        0.158472, 0.158123, 0.157804, 0.157485, 0.157166, 0.156847, 0.156528, 0.156209, 0.155890, 0.155571, 0.155252,
        0.154933, 0.154614, 0.154295, 0.153976, 0.153657, 0.153345, 0.153068, 0.152791, 0.152513, 0.152236, 0.151958,
        0.151681, 0.151403, 0.151097, 0.150774, 0.150450, 0.150127, 0.149803, 0.149479, 0.149156, 0.148835, 0.148548,
        0.148261, 0.147975, 0.147688, 0.147401, 0.147115, 0.146828, 0.146539, 0.146248, 0.145956, 0.145665, 0.145374,
        0.145083, 0.144791, 0.144500, 0.144276, 0.144054, 0.143832, 0.143610, 0.143388, 0.143166, 0.142944, 0.142717,
        0.142481, 0.142245, 0.142009, 0.141773, 0.141538, 0.141302, 0.141066, 0.140806, 0.140543, 0.140279, 0.140016,
        0.139752, 0.139489, 0.139225, 0.138966, 0.138716, 0.138467, 0.138217, 0.137967, 0.137718, 0.137468, 0.137218,
        0.137032, 0.136865, 0.136699, 0.136533, 0.136366, 0.136200, 0.136033, 0.135849, 0.135600, 0.135350, 0.135100,
        0.134851, 0.134601, 0.134352, 0.134102, 0.133907, 0.133740, 0.133574, 0.133408, 0.133241, 0.133075, 0.132908,
        0.132734, 0.132499, 0.132263, 0.132027, 0.131791, 0.131556, 0.131320, 0.131084, 0.130858, 0.130641, 0.130424,
        0.130206, 0.129989, 0.129772, 0.129555, 0.129338, 0.129236, 0.129134, 0.129032, 0.128931, 0.128829, 0.128727,
        0.128626, 0.128476, 0.128268, 0.128060, 0.127852, 0.127644, 0.127436, 0.127228, 0.127020, 0.126799, 0.126577,
        0.126355, 0.126133, 0.125911, 0.125689, 0.125468, 0.125270, 0.125117, 0.124965, 0.124812, 0.124659, 0.124507,
        0.124354, 0.124202, 0.124093, 0.123996, 0.123899, 0.123802, 0.123705, 0.123608, 0.123511, 0.123387, 0.123178,
        0.122970, 0.122762, 0.122554, 0.122346, 0.122138, 0.121930, 0.121761, 0.121608, 0.121455, 0.121303, 0.121150,
        0.120998, 0.120845, 0.120698, 0.120588, 0.120477, 0.120366, 0.120255, 0.120144, 0.120033, 0.119922, 0.119816,
        0.119715, 0.119613, 0.119511, 0.119409, 0.119308, 0.119206, 0.119101, 0.118926, 0.118750, 0.118574, 0.118399,
        0.118223, 0.118047, 0.117872, 0.117748, 0.117678, 0.117609, 0.117540, 0.117470, 0.117401, 0.117332, 0.117262,
        0.117102, 0.116936, 0.116769, 0.116603, 0.116436, 0.116270, 0.116103, 0.115948, 0.115809, 0.115670, 0.115532,
        0.115393, 0.115254, 0.115115, 0.114977, 0.114930, 0.114903, 0.114875, 0.114847, 0.114819, 0.114792, 0.114764,
        0.114712, 0.114596, 0.114480, 0.114365, 0.114249, 0.114134, 0.114018, 0.113902, 0.113733, 0.113543, 0.113354,
        0.113164, 0.112975, 0.112785, 0.112596, 0.112435, 0.112407, 0.112379, 0.112352, 0.112324, 0.112296, 0.112268,
        0.112241, 0.112109, 0.111914, 0.111720, 0.111526, 0.111332, 0.111138, 0.110944, 0.110764, 0.110759, 0.110754,
        0.110750, 0.110745, 0.110740, 0.110736, 0.110731, 0.110717, 0.110694, 0.110671, 0.110648, 0.110624, 0.110601,
        0.110578, 0.110555, 0.110348, 0.110135, 0.109922, 0.109710, 0.109497, 0.109284, 0.109072, 0.108945, 0.108935,
        0.108926, 0.108917, 0.108908, 0.108898, 0.108889, 0.108880, 0.108854, 0.108827, 0.108799, 0.108771, 0.108744,
        0.108716, 0.108688, 0.108608, 0.108413, 0.108219, 0.108025, 0.107831, 0.107637, 0.107442, 0.107248, 0.107203,
        0.107203, 0.107203, 0.107203, 0.107203, 0.107203, 0.107203, 0.107185, 0.107102, 0.107018, 0.106935, 0.106852,
        0.106769, 0.106686, 0.106602, 0.106482, 0.106344, 0.106205, 0.106066, 0.105928, 0.105789, 0.105650, 0.105527,
        0.105527, 0.105527, 0.105527, 0.105527, 0.105527, 0.105527, 0.105527, 0.105504, 0.105462, 0.105420, 0.105379,
        0.105337, 0.105296, 0.105254, 0.105212, 0.105073, 0.104934, 0.104795, 0.104657, 0.104518, 0.104379, 0.104241,
        0.104146, 0.104105, 0.104063, 0.104021, 0.103980, 0.103938, 0.103897, 0.103855, 0.103851, 0.103851, 0.103851,
        0.103851, 0.103851, 0.103851, 0.103851, 0.103841, 0.103814, 0.103786, 0.103758, 0.103730, 0.103703, 0.103675,
        0.103647, 0.103575, 0.103492, 0.103408, 0.103325, 0.103242, 0.103159, 0.103076, 0.102996, 0.102926, 0.102857,
        0.102788, 0.102718, 0.102649, 0.102580, 0.102510, 0.102460, 0.102419, 0.102377, 0.102336, 0.102294, 0.102252,
        0.102211, 0.102171, 0.102143, 0.102116, 0.102088, 0.102060, 0.102032, 0.102005, 0.101977, 0.101850, 0.101656,
        0.101462, 0.101268, 0.101073, 0.100879, 0.100685, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499,
        0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499,
        0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499, 0.100499,
        0.100499, 0.100499, 0.100499, 0.100503, 0.100508, 0.100512, 0.100517, 0.100522, 0.100526, 0.100531, 0.100481,
        0.100296, 0.100111, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000,
        0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000,
        0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000,
        0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000,
        0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000
    };
    return humanCai[index];
} // CBTensionModelLand17::CaBrixius

TFloat CBTensionModelLand17::CaCoppini(TFloat t) {
    // Calcium transient that was used in the Land17 publication.
    // Values taken from the Matlab script downloaded at http://cemrg.co.uk/models.html -> A MODEL OF CARDIAC CONTRACTION BASED ON NOVEL MEASUREMENTS OF TENSION DEVELOPMENT IN HUMAN CARDIOMYOCYTES.
    t = (t - GetActivationTime());
    while (t > cycleLength_) {
        t = t - cycleLength_;
    }
    int index = t/1e-3;
    if (index < 0) {
        index = 0;
    }
    
    // if cycleLength > 1s, use last index
    index = std::min(1000, index);
    double humanCai[] = {
        0.166, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.166, 0.170, 0.175, 0.181, 0.188, 0.196, 0.205, 0.215, 0.225,
        0.235, 0.247, 0.259, 0.270, 0.280, 0.289, 0.297, 0.304, 0.312, 0.319, 0.328, 0.336, 0.345, 0.353, 0.359, 0.364,
        0.368, 0.373, 0.377, 0.382, 0.388, 0.394, 0.401, 0.408, 0.413, 0.416, 0.420, 0.423, 0.427, 0.432, 0.438, 0.444,
        0.450, 0.455, 0.460, 0.463, 0.465, 0.468, 0.470, 0.475, 0.480, 0.485, 0.491, 0.496, 0.499, 0.502, 0.504, 0.504,
        0.506, 0.510, 0.513, 0.519, 0.524, 0.529, 0.531, 0.532, 0.532, 0.532, 0.533, 0.535, 0.538, 0.543, 0.547, 0.551,
        0.554, 0.554, 0.554, 0.554, 0.553, 0.555, 0.559, 0.561, 0.565, 0.568, 0.570, 0.569, 0.568, 0.568, 0.566, 0.567,
        0.570, 0.574, 0.577, 0.581, 0.581, 0.581, 0.579, 0.578, 0.578, 0.578, 0.579, 0.582, 0.586, 0.588, 0.588, 0.587,
        0.585, 0.583, 0.581, 0.580, 0.583, 0.587, 0.590, 0.593, 0.593, 0.592, 0.590, 0.586, 0.584, 0.585, 0.586, 0.590,
        0.592, 0.594, 0.593, 0.593, 0.590, 0.588, 0.586, 0.586, 0.587, 0.589, 0.592, 0.594, 0.593, 0.591, 0.589, 0.586,
        0.584, 0.583, 0.585, 0.587, 0.589, 0.591, 0.589, 0.588, 0.585, 0.582, 0.580, 0.579, 0.580, 0.582, 0.585, 0.587,
        0.586, 0.584, 0.581, 0.579, 0.578, 0.577, 0.578, 0.580, 0.581, 0.583, 0.581, 0.579, 0.576, 0.572, 0.570, 0.569,
        0.572, 0.574, 0.576, 0.578, 0.577, 0.574, 0.572, 0.569, 0.567, 0.566, 0.567, 0.569, 0.570, 0.570, 0.569, 0.567,
        0.563, 0.560, 0.558, 0.557, 0.558, 0.559, 0.561, 0.562, 0.563, 0.560, 0.557, 0.554, 0.551, 0.550, 0.550, 0.552,
        0.553, 0.555, 0.553, 0.550, 0.546, 0.544, 0.541, 0.541, 0.542, 0.543, 0.545, 0.545, 0.543, 0.541, 0.537, 0.534,
        0.531, 0.531, 0.530, 0.532, 0.534, 0.535, 0.533, 0.530, 0.527, 0.523, 0.520, 0.518, 0.518, 0.519, 0.521, 0.521,
        0.520, 0.518, 0.515, 0.512, 0.509, 0.508, 0.508, 0.508, 0.510, 0.510, 0.508, 0.505, 0.502, 0.499, 0.497, 0.496,
        0.496, 0.496, 0.496, 0.496, 0.496, 0.493, 0.490, 0.486, 0.484, 0.483, 0.483, 0.483, 0.484, 0.483, 0.481, 0.479,
        0.476, 0.472, 0.469, 0.467, 0.467, 0.467, 0.468, 0.469, 0.466, 0.464, 0.460, 0.457, 0.454, 0.452, 0.453, 0.453,
        0.454, 0.454, 0.453, 0.450, 0.447, 0.444, 0.441, 0.439, 0.439, 0.440, 0.441, 0.441, 0.440, 0.437, 0.433, 0.430,
        0.427, 0.426, 0.426, 0.426, 0.426, 0.426, 0.424, 0.421, 0.418, 0.415, 0.413, 0.411, 0.410, 0.411, 0.412, 0.411,
        0.409, 0.407, 0.404, 0.401, 0.399, 0.398, 0.397, 0.397, 0.397, 0.397, 0.396, 0.393, 0.390, 0.387, 0.384, 0.383,
        0.383, 0.383, 0.384, 0.384, 0.382, 0.380, 0.377, 0.374, 0.372, 0.370, 0.369, 0.370, 0.370, 0.370, 0.368, 0.366,
        0.364, 0.361, 0.358, 0.357, 0.357, 0.357, 0.357, 0.356, 0.355, 0.353, 0.350, 0.348, 0.345, 0.344, 0.343, 0.344,
        0.344, 0.344, 0.343, 0.341, 0.338, 0.336, 0.333, 0.332, 0.332, 0.332, 0.332, 0.332, 0.332, 0.330, 0.327, 0.325,
        0.323, 0.322, 0.321, 0.321, 0.321, 0.321, 0.320, 0.318, 0.316, 0.313, 0.311, 0.310, 0.310, 0.310, 0.311, 0.311,
        0.310, 0.308, 0.306, 0.303, 0.301, 0.301, 0.300, 0.300, 0.301, 0.301, 0.300, 0.298, 0.296, 0.294, 0.292, 0.291,
        0.291, 0.291, 0.292, 0.292, 0.291, 0.289, 0.287, 0.285, 0.284, 0.283, 0.282, 0.283, 0.283, 0.282, 0.281, 0.280,
        0.279, 0.277, 0.275, 0.274, 0.274, 0.274, 0.275, 0.274, 0.273, 0.272, 0.270, 0.269, 0.268, 0.267, 0.267, 0.267,
        0.267, 0.267, 0.266, 0.265, 0.263, 0.261, 0.260, 0.260, 0.260, 0.260, 0.260, 0.260, 0.259, 0.258, 0.257, 0.255,
        0.253, 0.253, 0.253, 0.253, 0.254, 0.254, 0.253, 0.252, 0.250, 0.248, 0.247, 0.246, 0.246, 0.247, 0.247, 0.247,
        0.247, 0.246, 0.244, 0.243, 0.241, 0.241, 0.241, 0.242, 0.242, 0.242, 0.242, 0.240, 0.239, 0.237, 0.236, 0.236,
        0.236, 0.236, 0.237, 0.237, 0.236, 0.236, 0.234, 0.233, 0.232, 0.231, 0.232, 0.232, 0.233, 0.232, 0.232, 0.231,
        0.229, 0.228, 0.227, 0.226, 0.226, 0.227, 0.228, 0.228, 0.227, 0.226, 0.225, 0.224, 0.223, 0.222, 0.222, 0.223,
        0.223, 0.223, 0.223, 0.222, 0.221, 0.220, 0.219, 0.218, 0.219, 0.219, 0.220, 0.220, 0.220, 0.219, 0.218, 0.217,
        0.216, 0.215, 0.215, 0.216, 0.216, 0.216, 0.216, 0.215, 0.214, 0.213, 0.212, 0.211, 0.212, 0.212, 0.213, 0.213,
        0.213, 0.212, 0.211, 0.209, 0.209, 0.209, 0.209, 0.209, 0.210, 0.210, 0.209, 0.208, 0.207, 0.206, 0.205, 0.205,
        0.205, 0.206, 0.206, 0.206, 0.206, 0.206, 0.205, 0.204, 0.203, 0.203, 0.203, 0.204, 0.204, 0.204, 0.204, 0.203,
        0.202, 0.201, 0.200, 0.200, 0.200, 0.201, 0.201, 0.201, 0.201, 0.200, 0.199, 0.199, 0.198, 0.197, 0.198, 0.198,
        0.199, 0.199, 0.199, 0.198, 0.197, 0.196, 0.195, 0.195, 0.196, 0.196, 0.197, 0.197, 0.197, 0.196, 0.195, 0.194,
        0.194, 0.193, 0.194, 0.194, 0.195, 0.195, 0.195, 0.194, 0.193, 0.192, 0.191, 0.191, 0.191, 0.192, 0.192, 0.193,
        0.192, 0.192, 0.191, 0.190, 0.190, 0.189, 0.190, 0.190, 0.191, 0.191, 0.191, 0.190, 0.189, 0.189, 0.188, 0.188,
        0.188, 0.188, 0.189, 0.189, 0.189, 0.188, 0.187, 0.187, 0.186, 0.186, 0.186, 0.187, 0.187, 0.188, 0.187, 0.187,
        0.186, 0.185, 0.185, 0.184, 0.185, 0.185, 0.186, 0.186, 0.186, 0.185, 0.184, 0.184, 0.183, 0.183, 0.183, 0.184,
        0.185, 0.185, 0.185, 0.184, 0.183, 0.182, 0.182, 0.181, 0.182, 0.182, 0.183, 0.184, 0.184, 0.183, 0.182, 0.181,
        0.181, 0.181, 0.181, 0.181, 0.182, 0.182, 0.182, 0.182, 0.181, 0.180, 0.179, 0.179, 0.180, 0.180, 0.181, 0.181,
        0.181, 0.180, 0.180, 0.179, 0.178, 0.178, 0.178, 0.179, 0.180, 0.180, 0.180, 0.179, 0.179, 0.178, 0.177, 0.177,
        0.177, 0.178, 0.179, 0.179, 0.179, 0.178, 0.178, 0.177, 0.176, 0.176, 0.176, 0.177, 0.178, 0.178, 0.178, 0.177,
        0.176, 0.175, 0.175, 0.175, 0.175, 0.176, 0.176, 0.177, 0.177, 0.177, 0.176, 0.175, 0.175, 0.175, 0.175, 0.176,
        0.176, 0.177, 0.176, 0.176, 0.175, 0.174, 0.174, 0.174, 0.174, 0.174, 0.175, 0.175, 0.175, 0.175, 0.174, 0.173,
        0.173, 0.173, 0.173, 0.174, 0.174, 0.174, 0.174, 0.174, 0.173, 0.172, 0.172, 0.172, 0.172, 0.172, 0.173, 0.173,
        0.173, 0.173, 0.172, 0.171, 0.171, 0.171, 0.171, 0.172, 0.172, 0.173, 0.173, 0.172, 0.171, 0.171, 0.170, 0.170,
        0.170, 0.171, 0.172, 0.172, 0.172, 0.172, 0.171, 0.170, 0.170, 0.169, 0.170, 0.171, 0.171, 0.171, 0.171, 0.171,
        0.170, 0.169, 0.169, 0.169, 0.169, 0.170, 0.170, 0.171, 0.171, 0.170, 0.170, 0.169, 0.168, 0.168, 0.168, 0.169,
        0.170, 0.170, 0.170, 0.170, 0.169, 0.168, 0.168, 0.167, 0.168, 0.169, 0.169, 0.170, 0.170, 0.169, 0.168, 0.168,
        0.167, 0.167, 0.168, 0.168, 0.169, 0.169, 0.169, 0.169, 0.168, 0.167, 0.167, 0.167, 0.167, 0.168, 0.168, 0.169,
        0.169, 0.168, 0.167, 0.167, 0.166, 0.167, 0.167, 0.168, 0.168, 0.168, 0.168, 0.168, 0.167, 0.166, 0.166, 0.166,
        0.166, 0.166, 0.167, 0.168, 0.167, 0.167, 0.166, 0.166, 0.165, 0.165, 0.166, 0.167, 0.167, 0.168, 0.168, 0.167,
        0.166, 0.166, 0.165, 0.165, 0.166, 0.166, 0.167, 0.167, 0.167, 0.167, 0.166, 0.165, 0.165, 0.165, 0.166, 0.166,
        0.167, 0.167, 0.167, 0.166, 0.166, 0.165, 0.165, 0.165, 0.165, 0.166, 0.166, 0.167, 0.167, 0.166, 0.166, 0.165,
        0.165, 0.165, 0.165, 0.166, 0.166, 0.167, 0.167, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.166, 0.166, 0.166,
        0.166, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.165, 0.166, 0.167, 0.167, 0.166, 0.165, 0.165, 0.164, 0.164,
        0.165, 0.165, 0.166, 0.166, 0.166, 0.166, 0.165, 0.165, 0.164, 0.164, 0.165, 0.166, 0.166, 0.166, 0.166, 0.166,
        0.165, 0.165, 0.164, 0.164, 0.164, 0.165, 0.166, 0.166, 0.166
    };
    return humanCai[index];
} // CBTensionModelLand17::CaCoppini

TFloat CBTensionModelLand17::CaExternal(TFloat t) {
    // Method to process a calcium transient from an external file
    t = (t - GetActivationTime());
    while (t > cycleLength_) {
        t = t - cycleLength_;
    }
    int index = t/1e-3;
    if (index < 0) {
        index = 0;
    }
    return humanCai_[index];
}

TFloat CBTensionModelLand17::CaElphy(TFloat t) {
    // Method to process calcium from electrophysiology
    TFloat calcium = 0;
    if (e_->GetType() == "T4")
        calcium = caiQP_[0];
    // even though we have calcium at QPs, there is currently only one tension model initialized per element
    // therefore, we use the centroid value as well in case of T10
    else if (e_->GetType() == "T10")
        calcium = caiQP_[0];
    return calcium;
}

TFloat CBTensionModelLand17::Cai(std::string flag, TFloat t) {
    if (flag == "Coppini")
        return CaCoppini(t);
    else if (flag == "Brixius")
        return CaBrixius(t);
    else if (flag == "External")
        return CaExternal(t);
    else if (flag == "Elphy")
        return CaElphy(t);
    else
        throw std::runtime_error("CBTensionModelLand17::Cai: unknown calcium transient option -> " + flag + ".");
}

TFloat CBTensionModelLand17::Overlap(TFloat lambda) {
    return 1.0 + beta0_ * (lambda + std::min(0.87, lambda) - 1.87);
}

/// Read a value from xml with fallback. Priorities: key+tail > keyFallback+tail > defaultValue
static TFloat InitKey(ParameterMap *params, std::string key, std::string keyFallback, std::string tail,
                      TFloat defaultValue) {
    if ((params->IsAvailable(key + tail) == false) && (params->IsAvailable(keyFallback + tail) == true))
        return params->Get<TFloat>(keyFallback);
    else
        return params->Get<TFloat>(key+tail, defaultValue);
}

CBTensionModelLand17::CBTensionModelLand17(CBElementSolid *e, ParameterMap *params) {
    // set model parameters
    Matrix3<TFloat> deformationTensor;
    e_ = e;
    assert(e_ != nullptr);
    caiQP_.resize(e_->GetNumberOfQuadraturePoints());
    e->GetDeformationTensor(deformationTensor);
    Tmax_ = e->GetMaterial()->GetProperties()->tensionMax_;
    int mi = e->GetMaterialIndex();
    ei_ = e->GetIndex();
    std::string key = "Materials.Mat_" + std::to_string(mi);
    std::string keyFallback = "Materials.Mat_Default";
    InitParamsFromXml(params, key+".Land17", keyFallback+".Land17");
    
    // set initial values for the state variables
    S_.t          = 0.0;
    S_.delta_t    = 0.0;
    S_.XS         = 0.0;
    S_.XW         = 0.0;
    S_.TRPN       = 0.0;
    S_.TmBlocked  = 1.0;
    S_.ZetaS      = 0.0;
    S_.ZetaW      = 0.0;
    S_.Cd         = 0.0;
    S_.lambda     = sqrt(deformationTensor.GetCol(0)*deformationTensor.GetCol(0));
    S_.dlambdadt  = 0.0;
    S_.Ta         = 0.0;
    S_.Tension    = 0.0;
    
    S_curr_ = S_;
    S_prev_ = S_;
    
    // set constant values that are dependent on init parameters
    kws_ = xi_ * kws_;
    cw_ = phi_ * kuw_ * ((1.0 - rs_) * (1.0 - rw_)) / ((1.0 - rs_) * rw_);
    cs_ = phi_ * kws_ * ((1.0 - rs_) * rw_) / rs_;
    kwu_ = kuw_ * (1.0 / rw_ - 1.0) - kws_;
    ksu_ = (kws_ * rw_ * (1.0 / rs_ - 1.0));
    A_ = Aeff_ * rs_ / ((1.0 - rs_) * rw_ + rs_);
    XSSS_ = rs_ * 0.5;
    XWSS_ = (1.0 - rs_) * rw_ * 0.5;
    ktm_block_ = ku_ * pow(TRPN50_, nTM_) * 0.5 / (0.5 - XSSS_ - XWSS_);
}

void CBTensionModelLand17::SaveStateVariablesAsNeeded(TFloat time) {
    // time is bigger than in the last call? -> the last call was successful, store S_curr to S_prev
    // time is smaller than in last call? -> previous run was non-successful, discard last computed values S_curr
    // time is roughly the same? -> update the current guess S_curr (with what?) ...
    if (time - S_curr_.t > 0) {
        S_prev_ = S_curr_;
        
        // export state variables on a successful export
        if (doExport_)
            this->WriteToFile(S_prev_);
    } else if (time - S_prev_.t < 0) {
        S_curr_ = S_prev_;
    } else {
        S_curr_ = S_prev_;
    }
    S_ = S_prev_;
}

double CBTensionModelLand17::CalcActiveTension(const math_pack::Matrix3<double> &deformation, const double time) {
    SaveStateVariablesAsNeeded(time);
    S_.t = time; // in s
    S_.delta_t = (S_.t - S_prev_.t) * 1000.0; // in s -> ms
    
    // compute / integrate variables
    // update extention ratio (lambda) and exention rate (dlambdadt)
    S_.lambda = sqrt(deformation.GetCol(0)*deformation.GetCol(0)); // in m/m
    if ((S_.delta_t == 0.0) || (rateDependancy_ == "OFF")) { // Avoid division by 0 in the first time-step
        S_.dlambdadt = 0.0;
    } else {
        S_.dlambdadt = ((S_.lambda - S_prev_.lambda) / S_.delta_t); // 1/ms
    }
    
    // XB model
    TFloat lambda_m = std::min(1.2, S_.lambda);
    TFloat h = std::max(0.0, Overlap(lambda_m));
    
    // unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
    TFloat XU = (1.0 - S_prev_.TmBlocked) - S_prev_.XW - S_prev_.XS;
    TFloat xb_ws = kws_ * S_prev_.XW;
    TFloat xb_uw = kuw_ * XU;
    TFloat xb_wu = kwu_ * S_prev_.XW;
    TFloat xb_su = ksu_ * S_prev_.XS;
    
    TFloat zs_pos = (S_prev_.ZetaS > 0.0) * S_prev_.ZetaS;
    TFloat zs_neg = (S_prev_.ZetaS < -1.0) * (-S_prev_.ZetaS-1.0);
    TFloat zs = std::max(zs_pos, zs_neg);
    TFloat gamma_rate = gammas_ * zs;
    TFloat xb_su_gamma = gamma_rate * S_prev_.XS;
    S_.XS = S_prev_.XS + S_.delta_t * dXSdt(xb_ws, xb_su, xb_su_gamma);
    
    TFloat gr_w = ((S_prev_.ZetaW < 0.0) ? -S_prev_.ZetaW : S_prev_.ZetaW); // use absolute value of ZetaW
    TFloat gamma_rate_w = gammaw_ * gr_w; // weak xbs don't like being strained
    TFloat xb_wu_gamma = gamma_rate_w * S_prev_.XW;
    S_.XW = S_prev_.XW + S_.delta_t * dXWdt(xb_uw, xb_wu, xb_ws, xb_wu_gamma);
    
    TFloat Ca50 = Ca50_ + beta1_ * (lambda_m - 1.0);
    S_.TRPN = S_prev_.TRPN + S_.delta_t * dTRPNdt(TRPNk_, Cai(calciumTransientType_, S_.t), Ca50, TRPNn_, S_prev_.TRPN);
    
    TFloat TRPN_NP = pow(S_.TRPN, (-nTM_/2.0));
    TFloat TRPN_NP_ = std::min(100.0, TRPN_NP);
    S_.TmBlocked = S_prev_.TmBlocked + S_.delta_t * dTmBlockeddt(ktm_block_, TRPN_NP_, XU, ku_, S_.TRPN, nTM_,
                                                                 S_prev_.TmBlocked);
    
    // velocity dependance - assumes distortion resets on W->S
    S_.ZetaS = S_prev_.ZetaS + S_.delta_t * dZetaSdt(A_, S_.dlambdadt, cs_, S_prev_.ZetaS);
    S_.ZetaW = S_prev_.ZetaW + S_.delta_t * dZetaWdt(A_, S_.dlambdadt, cw_, S_prev_.ZetaW);
    
    // ActiveTension in kPa
    // Tmax_ = 1000 should be used to do kPa -> Pa
    // don't allow negative values
    S_.Ta = std::max(0.0, h * (Tref_ / rs_) * ((S_.ZetaS + 1.0) * S_.XS + S_.ZetaW * S_.XW));
    
    // Minimal implementation of the passive cell model
    // Similar to a standard linear solid model. It is used for the viscoelastic response.
    TFloat C_s = (S_.lambda - 1.0) - S_prev_.Cd;
    TFloat eta = (C_s > 0.0) ? eta_l_ : eta_s_;
    S_.Cd = S_prev_.Cd + S_.delta_t * dC_ddt(k_, eta, C_s);
    TFloat F_d = a_ * k_ * C_s;
    
    // Total Tension
    // the sum should be non-negative as well
    S_.Tension = std::max(0.0, S_.Ta + F_d);
    
    S_curr_ = S_;
    
    return Tmax_ * S_.Tension;
    
    //  return Tmax_ * (S_.Tension + CalcActiveStiffness() * (S_.lambda - S_prev_.lambda));
} // CBTensionModelLand17::CalcActiveTension

void CBTensionModelLand17::InitParamsFromXml(ParameterMap *parameters, std::string parameterKey,
                                             std::string parameterKeyFallback) {
    rateDependancy_ = parameters->Get<std::string>(parameterKey + ".rateDependancy", "OFF");
    startTime_    = InitKey(parameters, parameterKey, parameterKeyFallback, ".startTime", 0.0);
    cycleLength_  = InitKey(parameters, parameterKey, parameterKeyFallback, ".cycleLength", 1.0);
    kff_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".kff", 1.0);
    kss_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".kss", 0.0);
    knn_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".knn", 0.0);
    ksn_ = InitKey(parameters, parameterKey, parameterKeyFallback, ".ksn", 0.0);
    
    SetStressCoefficients(Matrix3<TFloat> {kff_, 0, 0,  0, kss_, ksn_,  0, ksn_, knn_});
    AddToActivationTime(startTime_);
    
    calciumTransientType_ = parameters->Get<std::string>(parameterKey + ".CalciumTransientType", "Coppini");
    calciumFile_          = parameters->Get<std::string>(parameterKey + ".CalciumFile", "");
    if (!(calciumTransientType_ == "Coppini") && !(calciumTransientType_ == "Brixius") && !(calciumTransientType_ == "External") && !(calciumTransientType_ == "Elphy")) {
        throw std::runtime_error(
                                 "CBTensionModelLand17::InitParamsFromXml: no suitable calcium transient was defined by the user (Coppini, Brixius, External (from file), Elphy (coupled to electrophysiology); default NONE");
    }
    if (calciumTransientType_ == "External") {
        ReadExternalCai(calciumFile_);
    }
    
    TRPNn_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".TRPNn", 2.0);
    TRPNk_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".TRPNk", 0.1);
    Ca50_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".Ca50", 0.805);
    ku_           = InitKey(parameters, parameterKey, parameterKeyFallback, ".ku", 1.0);
    nTM_          = InitKey(parameters, parameterKey, parameterKeyFallback, ".nTM", 5.0);
    TRPN50_       = InitKey(parameters, parameterKey, parameterKeyFallback, ".TRPN50", 0.35);
    kuw_          = InitKey(parameters, parameterKey, parameterKeyFallback, ".kuw", 0.182);
    kws_          = InitKey(parameters, parameterKey, parameterKeyFallback, ".kws", 0.012);
    rs_           = InitKey(parameters, parameterKey, parameterKeyFallback, ".rs", 0.25);
    rw_           = InitKey(parameters, parameterKey, parameterKeyFallback, ".rw", 0.5);
    gammas_       = InitKey(parameters, parameterKey, parameterKeyFallback, ".gammas", 0.0085);
    gammaw_       = InitKey(parameters, parameterKey, parameterKeyFallback, ".gammaw", 0.615);
    phi_          = InitKey(parameters, parameterKey, parameterKeyFallback, ".phi", 2.23);
    Aeff_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".Aeff", 25.0);
    beta0_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".beta0", 2.3);
    beta1_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".beta1", -2.4);
    Tref_         = InitKey(parameters, parameterKey, parameterKeyFallback, ".Tref", 120.0);
    a_            = InitKey(parameters, parameterKey, parameterKeyFallback, ".a", 2.1);
    k_            = InitKey(parameters, parameterKey, parameterKeyFallback, ".k", 7.0);
    eta_l_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".eta_l", 200.0);
    eta_s_        = InitKey(parameters, parameterKey, parameterKeyFallback, ".eta_s", 20.0);
    xi_           = InitKey(parameters, parameterKey, parameterKeyFallback, ".xi", 1.0);
    
    // init variables for exporting
    std::vector<TInt> elementIndicesVec = parameters->GetArray<TInt>(parameterKey+".ExportIndices", std::vector<TInt>());
    std::set<TInt> elementIndicesSet(elementIndicesVec.begin(), elementIndicesVec.end());
    auto search = elementIndicesSet.find(ei_);
    if (search != elementIndicesSet.end())
        doExport_ = true;
    if ((parameters->IsAvailable(parameterKey + ".Filename") == false) &&
        (parameters->IsAvailable(parameterKeyFallback + ".Filename") == true) )
        filename_ = parameters->Get<std::string>(parameterKeyFallback);
    else
        filename_ = parameters->Get<std::string>(parameterKey+".Filename", "Land17.dat");
} // CBTensionModelLand17::InitParamsFromXml

inline void CBTensionModelLand17::WriteToFile(const StateVariables &S) {
    std::ofstream file;
    
    file.open(filename_.c_str(), std::ios::app);
    if (!file.good())
        throw std::runtime_error("CBTensionModelLand17::WriteToFile: Couldn't create " + filename_ + ".");
    
    // write header
    if (!headerWritten_) {
        file << "index";
        file << "\t" << "t";
        file << "\t" << "delta_t";
        file << "\t" << "XS";
        file << "\t" << "XW";
        file << "\t" << "TRPN";
        file << "\t" << "TmBlocked";
        file << "\t" << "ZetaS";
        file << "\t" << "ZetaW";
        file << "\t" << "lambda";
        file << "\t" << "dlambdadt";
        file << "\t" << "Cai";
        file << "\t" << "Ta";
        file << "\t" << "Overlap";
        file << "\t" << "Tension";
        file << "\t" << std::endl;
        headerWritten_ = true;
    }
    
    // write content to file
    file << std::setprecision(14); // needed to increase precision of export to compare with matlab
    file << ei_;
    file << "\t" << S.t;
    file << "\t" << S.delta_t;
    file << "\t" << S.XS;
    file << "\t" << S.XW;
    file << "\t" << S.TRPN;
    file << "\t" << S.TmBlocked;
    file << "\t" << S.ZetaS;
    file << "\t" << S.ZetaW;
    file << "\t" << S.lambda;
    file << "\t" << S.dlambdadt;
    file << "\t" << Cai(calciumTransientType_, S.t);
    file << "\t" << S.Ta;
    file << "\t" << Overlap(std::min(1.2, S.lambda));
    file << "\t" << S.Tension;
    file << std::endl;
    file.close();
} // CBTensionModelLand17::WriteToFile

inline void CBTensionModelLand17::ReadExternalCai(std::string filename) {
    std::ifstream file(filename);
    
    for (double a; file >> a;) {
        humanCai_.push_back(a);
    }
}

CBStatus CBTensionModelLand17::SetActiveTensionAtQuadraturePoint(TInt indexQP, TFloat Cai) {
    CBStatus rc = CBStatus::SUCCESS;
    
    if (e_->GetType() == "T4") {
        if (indexQP == 0) {
            caiQP_[indexQP] = Cai;
        } else {
            throw std::runtime_error("CBTensionModelLand17::SetActiveTensionAtQuadraturePoint: try to access a not existing quadrature point for T4 elements.");
        }
    } else if (e_->GetType() == "T10") {
        if (indexQP < 5) {
            caiQP_[indexQP] = Cai;
        } else {
            throw std::runtime_error("CBTensionModelLand17::SetActiveTensionAtQuadraturePoint: try to access a not existing quadrature point for T10 elements.");
        }
    } else {
        rc = CBStatus::FAILED;
    }
    
    return rc;
}
