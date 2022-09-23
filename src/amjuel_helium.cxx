#include "../include/amjuel_helium.hxx"

////////////////////////////////////////////////////////////
// e + he -> he+ + 2e

/// Coefficients to calculate the reaction rate <σv>
/// Amjuel reaction 2.3.9a, page 161
static constexpr const BoutReal he01_rate_coefs[9][9] = {
    {-42.27118452798, 0.1294554451998, -0.08433979538052, 0.04910721979375,
     -0.01454047282438, 0.002178105605879, -0.0001657512355348, 6.161429564793e-06,
     -8.910615590909e-08},
    {24.11668100975, -0.08121999208281, 0.04052570160482, -0.02367924962508,
     0.008488392041366, -0.001452752408581, 0.0001170902182939, -4.410479245308e-06,
     6.297315949647e-08},
    {-12.03181133667, -0.003998282970932, -0.00281991919306, -0.00190488772724,
     -0.0002390948585334, 0.0001844484422285, -1.97272802786e-05, 7.779440219801e-07,
     -1.033814145233e-08},
    {3.829444688521, 0.02546414073266, 0.002654490306111, 0.001087493205419,
     -0.0004469192206896, 3.71553815559e-05, -1.595144154431e-06, 6.311039124056e-08,
     -1.48598916668e-09},
    {-0.7945839257175, -0.0149359787485, -0.001018320076497, 0.0002821927325759,
     3.269264854581e-05, -5.937518354028e-06, 4.714656637197e-07, -2.433462923993e-08,
     5.307423532159e-10},
    {0.1054334178555, 0.004338821244147, -0.0001483560478208, -6.901574689672e-05,
     6.350490312899e-06, -4.414167358057e-07, 1.266603603049e-08, 8.049435558339e-10,
     -3.807796193572e-11},
    {-0.008578643565653, -0.0006689202603525, 9.084162487421e-05, -4.184111347149e-06,
     1.153919327151e-07, 3.797435455934e-08, -4.123383037275e-09, 1.095960078746e-10,
     -5.109801608123e-14},
    {0.0003886232727181, 5.180805123476e-05, -1.125453787291e-05, 1.536214841434e-06,
     -1.632601398517e-07, 8.948177075796e-09, -1.853674996294e-10, 1.342166707999e-14,
     1.184569645146e-14},
    {-7.487575233223e-06, -1.58297743374e-06, 4.413792107083e-07, -7.832095176637e-08,
     9.58697477495e-09, -6.73907617081e-10, 2.565598443992e-11, -4.994625098807e-13,
     4.12404880445e-15}};

/// Effective electron cooling rate due to ionization of Helium atoms.
/// Fujimoto Formulation II (only ground level transported, no meta-stables kept
/// explicitly)
static constexpr const BoutReal he01_radiation_coefs[9][9] = {
    {-35.35258393674, -0.03428249311738, 0.06378071832382, -0.02849818870377,
     0.006041903480645, -0.000686453216556, 4.251155616815e-05, -1.351759350582e-06,
     1.728801977101e-08},
    {19.81855871044, 0.04854482688892, -0.05088928946831, 0.01732110218818,
     -0.002781419068092, 0.0002244804771683, -8.875290574348e-06, 1.399429819761e-07,
     -1.38977874051e-10},
    {-9.334355651224, -0.04524206463148, 0.02103002869692, -0.004463941003028,
     0.0002900917070658, 2.482449118881e-05, -4.278064413224e-06, 2.040570181783e-07,
     -3.324224092217e-09},
    {2.80031425041, 0.0247435078798, -0.006012991773715, 0.0008918009845745,
     -2.616249899141e-05, -6.885545577757e-06, 7.013616309712e-07, -2.570063437935e-08,
     3.573487194914e-10},
    {-0.5489088598705, -0.007339538872774, 0.0007783071302508, -4.483274558979e-05,
     1.900991581685e-06, -9.747171692727e-07, 1.349829568374e-07, -5.815812094637e-09,
     6.686532777575e-11},
    {0.06902095610357, 0.001234159378604, 2.989745411104e-05, -3.04090620334e-05,
     2.951386149372e-06, 7.592185107575e-08, -1.805060230413e-08, 3.156859219121e-10,
     1.07116869734e-11},
    {-0.00534294006913, -0.0001223169549107, -1.500790305823e-05, 5.253922160283e-06,
     -4.468905893926e-07, 7.483496971361e-09, -9.777558713428e-10, 1.770619394125e-10,
     -6.050995244427e-12},
    {0.0002313175089975, 6.966436907981e-06, 8.94496290981e-07, -1.712024596447e-07,
     -9.782015167261e-09, 2.499416349949e-09, 4.731973382221e-11, -1.845161957843e-11,
     6.01107014323e-13},
    {-4.279800193256e-06, -1.81546666991e-07, -2.282174576618e-09, -6.972920569943e-09,
     2.60719149454e-09, -2.870919514967e-10, 8.059675146168e-12, 3.704316808942e-13,
     -1.713225271579e-14}};

void AmjuelHeIonisation01::transform(Options& state, Field3D &reaction_rate) {
  electron_reaction(state["species"]["e"],
                    state["species"]["he"],  // From helium atoms
                    state["species"]["he+"], // To helium ions
                    he01_rate_coefs, he01_radiation_coefs,
                    0.0, // Note: Ionisation potential included in radiation_coefs,
                    reaction_rate

  );
}

////////////////////////////////////////////////////////////
// e + he+ -> he

/// Coefficients to calculate the reaction rate <σv>
/// Amjuel reaction 2.3.13a
static constexpr const BoutReal he10_rate_coefs[9][9] = {
    {-28.72754373123, -0.006171082987797, 0.02414548639597, -0.007188662067622,
     0.0009481268604767, -1.958887458637e-05, -5.507786383328e-06, 4.35828868693e-07,
     -9.50327209101e-09},
    {1.564233603544, -0.03972220721457, -0.04466712599181, 0.01247359158796,
     -0.001660591942878, 6.019181402025e-05, 3.800156798817e-06, -3.377807793756e-07,
     6.828447501225e-09},
    {-6.182140631482, 0.1626641668186, 0.03366589582541, -0.007413737965595,
     0.001220189896183, -9.50529572475e-05, 4.459492214068e-06, -1.552772441333e-07,
     2.866586118879e-09},
    {5.459428677778, -0.1700323494998, -0.01540106384088, 0.0009524545793262,
     -8.734341535385e-05, -2.796027477899e-06, 4.561981097438e-07, 4.940311502014e-09,
     -6.52572501076e-10},
    {-2.128115924661, 0.07233939709414, 0.005819196258503, -7.655935845761e-05,
     -1.83794906705e-05, 4.72578983298e-06, -3.99778241186e-07, 1.036731541123e-08,
     -3.373845712183e-11},
    {0.4373730373037, -0.01574917019835, -0.001456253436544, 4.772491845078e-05,
     -1.827059132463e-06, -6.94116329271e-08, 2.716740135949e-08, -1.143121626264e-09,
     1.295139027087e-11},
    {-0.04972257208732, 0.001866175274689, 0.0002047337498511, -1.004438052808e-05,
     7.59073486585e-07, -6.771179147667e-08, 1.218720257518e-09, 6.78702447954e-11,
     -2.432253541918e-12},
    {0.002967287371427, -0.0001147811325052, -1.460813593905e-05, 7.422385993164e-07,
     -3.281946488134e-08, 2.164459880579e-09, 1.113868237282e-10, -1.513922678655e-11,
     3.951084520871e-13},
    {-7.271204747116e-05, 2.874049670122e-06, 4.124421172202e-07, -1.689203971933e-08,
     -9.071172814458e-10, 1.844295219334e-10, -2.055023511556e-11, 1.101902611511e-12,
     -2.206082129473e-14}};

/// Radiation energy loss from helium recombination
/// The potential energy (24.586eV per event) should be added to the electrons
/// so that the process may be a net energy source for the electrons
static constexpr const BoutReal he10_radiation_coefs[9][9] = {
    {-25.38377692766, -0.04826880987619, 0.0679657596731, -0.0240171002139,
     0.004130156138736, -0.0003494803122018, 1.34502510054e-05, -1.323917127568e-07,
     -2.551716207606e-09},
    {2.472758419513, 0.1668058989207, -0.1265192781981, 0.02938171777028,
     -0.00221652505507, -0.0001261523686946, 2.701401918133e-05, -1.370446267883e-06,
     2.313673787201e-08},
    {-8.864417999957, -0.1882326730037, 0.119402867431, -0.02382836629119,
     0.001134820469638, 0.0001782113978272, -2.305554399898e-05, 9.43029409318e-07,
     -1.305188423829e-08},
    {8.394970578944, 0.08397993216045, -0.0579697281374, 0.01158600348753,
     -0.0007504743150582, -2.602911694939e-05, 4.568209602293e-06, -1.458110560501e-07,
     9.826599911934e-10},
    {-3.465864794112, -0.0157268418022, 0.01398192327776, -0.002700181027443,
     0.0001902304157269, -1.534387905925e-06, -1.03206007926e-07, -1.355858638619e-08,
     5.917279771473e-10},
    {0.7479071085372, 0.0005997666028811, -0.001614053457119, 0.0002620866439317,
     -1.103039382799e-05, -4.215447554819e-07, -9.926133276192e-09, 4.148813674084e-09,
     -1.207867670158e-10},
    {-0.08863575102304, 0.0001901540166344, 6.941090299375e-05, -3.042043168371e-06,
     -1.677907209787e-06, 2.153652742395e-07, -6.354480058307e-09, -1.223103792568e-10,
     5.54305794673e-12},
    {0.005484926807853, -2.510359436743e-05, 1.123735445147e-06, -9.198494797723e-07,
     2.058851315121e-07, -1.866416375894e-08, 7.564835556537e-10, -1.622993472948e-11,
     2.428986170198e-13},
    {-0.0001388441945179, 9.1419955967e-07, -1.16891589033e-07, 3.485370731777e-08,
     -5.086412415216e-09, 3.674153797642e-10, -1.621809988343e-11, 6.737654534264e-13,
     -1.678705755876e-14}};

void AmjuelHeRecombination10::transform(Options& state, Field3D &reaction_rate) {
  electron_reaction(state["species"]["e"],
                    state["species"]["he+"], // From helium ions
                    state["species"]["he"],  // To helium atoms
                    he10_rate_coefs, he10_radiation_coefs,
                    24.586, // Ionisation potential heating of electrons
                    reaction_rate
  );
}
