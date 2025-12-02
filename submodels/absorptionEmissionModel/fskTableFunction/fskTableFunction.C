#include "fskTableFunction.H"
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

void loadIntoMem(std::vector<float>& kq, std::string fpath)
{
    // FSK Table Parameters
    const int NDB = 32, Nq=32;
    const int nTgDB=28, nxco2DB=13, nxh2oDB=13, nxcoDB=6, nfvDB=2;
    const double xCO2DB[nxco2DB]={0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    const double xH2ODB[nxh2oDB]={0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    const double xCODB[nxcoDB]={0.e0,0.01e0,0.05e0,0.1e0,0.25e0,0.5e0};
    const double fvDB[nfvDB]={0.e0,1e-5};
    int nT=nTgDB;
    // Needs to be type float for right data size
    //float* kq = new(nothrow) float[nfvDB*nxco2DB*nxh2oDB*nxcoDB*nTgDB*nT*Nq];
    kq.assign(nfvDB*nxco2DB*nxh2oDB*nxcoDB*nTgDB*nT*Nq,0);
    float kDB[NDB + 1];

    cout << "Loading FSK Table v3 data into memory." << endl;
    double P = 1.0;
    for (int ifv = 0; ifv < nfvDB; ++ifv) //fv Loop
    {
        for (int ixco2 = 0; ixco2 < nxco2DB; ++ixco2) //CO2 Loop
        {
            for (int ixh2o = 0; ixh2o < nxh2oDB; ++ixh2o) //H2O Loop
            {
                for (int ixco = 0; ixco < nxcoDB; ++ixco) //CO Loop
                {
                    if ((xCO2DB[ixco2]+xH2ODB[ixh2o]+xCODB[ixco])>1.e0)
                    {
                        if ((xCO2DB[std::max(ixco2-1,1)]+xH2ODB[std::max(ixh2o-1,1)]+xCODB[std::max(ixco-1,1)])>=1.e0)
                        {
                            continue;
                        }
                    }

                    char buffer[100];
                    std::snprintf(buffer, 100, "P.%04.1f_CO2.%4.2f_H2O.%4.2f_CO.%4.2f_fv.%05.E_k.dat", P, xCO2DB[ixco2], xH2ODB[ixh2o], xCODB[ixco], fvDB[ifv]);
                    std::string saveNamek(buffer);
                    std::string dot=".";
                    saveNamek.insert(37, dot);
                    std::string filepath = fpath + saveNamek;
                    fstream file;
                    file.open(filepath, ios::in | ios::binary);
                    if (!file)
                    {
                        cout << "Error opening file! " << filepath << endl;
                        break;
                    }
                    else
                    {
                        for (int iTg = 0; iTg < nTgDB; ++iTg)
                        {
                            int recl = 4*(NDB + 1);
                            int id = nTgDB*(iTg)*recl;
                            for (int iT = 0; iT < nT; ++iT)
                            {

                                file.seekg(id);
                                file.read(reinterpret_cast<char*>(&kDB), 4 * (NDB + 1));
                                for (int iq = 0; iq < NDB; iq++)
                                {
                                    int index = iq + Nq*(iT + nT*(iTg + nTgDB*(ixco + nxcoDB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(ifv)))))); // Working! ! !
                                    kq[index] = kDB[iq];
                                    //cout << "Index " << index << " value is " << kDB[iq] << endl;
                                }
                                id = id + recl;
                            }
                        }
                    }
                    file.close();
                    //cout << "File " << filepath << " succesfully read!" << endl;
                }
            }
        }
    }
    cout << "FSK Table v3 data successfully loaded into memory." << endl;
}

void loadIntoMemV2(std::vector<float>& kq, std::string fpath, const double Tref)
{
    // FSK Table Parameters
    const int NDB = 32, Nq=32;
    const int nTgDB=28, nxco2DB=13, nxh2oDB=13, nxcoDB=6, nxnh3DB = 13, nxch4DB = 13;
    const double xCO2DB[nxco2DB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xH2ODB[nxh2oDB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xNH3DB[nxnh3DB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xCH4DB[nxch4DB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xCODB[nxcoDB]={0.0e0,0.010e0,0.050e0,0.10e0,0.25e0,0.50e0};
    int nT=nTgDB;
    // Needs to be type float for right data size
    //float* kq = new(nothrow) float[nfvDB*nxco2DB*nxh2oDB*nxcoDB*nTgDB*nT*Nq];
    // Loading only one reference temperature in this DataBase due to more mixture states compared to Table v3
    //kq.assign(nxco2DB*nxh2oDB*nxnh3DB*nxch4DB*nxcoDB*nT*Nq,0);
    // Loading only FSK Table v3 gases into memory
    //kq.assign(nxco2DB*nxh2oDB*nxcoDB*nTgDB*nT*Nq,0);
    // Loading only fuels
    kq.assign(nxch4DB*nxnh3DB*nTgDB*nT*Nq, 0);
    //float kDB[NDB + 1];

    cout << "Loading FSK Table (Ammonia DataBase) data into memory." << endl;
    double P = 1.0;

    //int index = 31 + Nq*(5 + nxcoDB*(12 + nxnh3DB*(12 + nxch4DB*(12 + nxh2oDB*(12 + nxco2DB*(27)))))); // Working! ! !
    //cout << "Largest index for Table: " << index << endl;
    for (int iTref = 0; iTref < nTgDB; ++iTref) // Tref Loop
    {
        double Tref_loop = 300.0 + iTref * 100.0;
        int TrefInt = static_cast<int>(std::round(Tref_loop));
        //int TrefInt = static_cast<int>(std::round(Tref));
        std::string filepath = fpath + "fsck_Tref_" + std::to_string(TrefInt) + ".bin";
        //cout << "Loading file: " << filepath << endl;
        fstream file;
        file.open(filepath, ios::in | ios::binary);
        if (!file)
        {
            cout << "Error opening file! " << filepath << endl;
            exit(EXIT_FAILURE);
        }
        // Open zero concentration state file
        std::string zeroFilePath = fpath + "zero_concentration_state_" + std::to_string(TrefInt) + ".bin";
        //cout << "Zero concentration state file : " << zeroFilePath << endl;
        std::ifstream zeroFile(fpath + "zero_concentration_state_" + std::to_string(TrefInt) + ".bin", std::ios::binary);
        if (!zeroFile)
        {
            cout << "Error opening zero concentration state file!" << endl;
            exit(EXIT_FAILURE);
        }
        int mixtureStates = 0;
        int overflowStates = 0;
        //int userCO2;
        //cout << "Enter the number of CO2 states to load (1-13): ";
        //cin >> userCO2;
        for (int iT = 0; iT < nT; ++iT) //fv Loop
        {
            for (int ixco2 = 0; ixco2 < 1; ++ixco2) //CO2 Loop
            {
                for (int ixh2o = 0; ixh2o < 1; ++ixh2o) //H2O Loop
                {
                    for (int ixch4 = 0; ixch4 < nxch4DB; ++ixch4) //CH4 Loop
                    {
                        for (int ixnh3 = 0; ixnh3 < nxnh3DB; ++ixnh3) //NH3 Loop
                        {
                            for (int ixco = 0; ixco < 1; ++ixco) //CO Loop
                            {
                                double xSum = xCO2DB[ixco2] + xH2ODB[ixh2o] + xNH3DB[ixnh3] + xCH4DB[ixch4] + xCODB[ixco];
                                if (xSum > 1.0e0)
                                {
                                    overflowStates++;
                                    //cout << ixco2 << "," << ixh2o << "," << ixch4 << "," << ixnh3 << "," << ixco << endl;
                                    continue;
                                }
                                else if (xSum < 1e-10)
                                {
                                    std::vector<float> k_zero(Nq, 0.0f);
                                
                                    zeroFile.read(reinterpret_cast<char*>(k_zero.data()), Nq * sizeof(float));
                                    for (int iq = 0; iq < Nq; ++iq)
                                    {
                                        // Only load CO, CO2, H2O gases
                                        //int index = iq + Nq*(ixco + nxcoDB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(iT + nT*(iTref)))));
                                        // Only load CH4 and NH3 gases
                                        int index = iq + Nq*(ixnh3 + nxnh3DB*(ixch4 + nxch4DB*(iT + nT*(iTref))));

                                        // Load all gases (8gb+)
                                        //int index = iq + Nq*(ixco + nxcoDB*(ixnh3 + nxnh3DB*(ixch4 + nxch4DB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(iT)))))); // Working! ! !
                                        
                                        kq[index] = k_zero[iq];
                                    }                               
                                    continue;
                                }

                                // Weird MATLAB precision issue requires skipping these two states (NEEDS FIX)
                                /*if ((ixco2 == 7 && ixh2o == 8 && ixch4 == 10 && ixnh3 == 5 && ixco == 3) || (ixco2 == 8 && ixh2o == 7 && ixch4 == 10 && ixnh3 == 5 && ixco == 3))
                                {
                                    overflowStates++;
                                    continue;
                                }*/

                                for (int iq = 0; iq < NDB; iq++)
                                {
                                    float value;
                                    file.read(reinterpret_cast<char*>(&value), sizeof(float));
                                    
                                    // Only load CO, CO2, H2O gases
                                    //if (ixch4>0 || ixnh3>0)
                                    // Only load CH4 and NH3 gases
                                    if (ixco > 0 || ixco2 > 0 || ixh2o > 0)
                                    {
                                        continue; // Skip CH4 and NH3 for FSK Table v3
                                    }
                                    // Only load CO, CO2, H2O gases
                                    //int index = iq + Nq*(ixco + nxcoDB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(iT + nT*(iTref)))));
                                    // Only load CH4 and NH3 gases
                                    int index = iq + Nq*(ixnh3 + nxnh3DB*(ixch4 + nxch4DB*(iT + nT*(iTref))));
                                    kq[index] = value;
                                    
                                    // Load all gases (8gb+)
                                    //int index = iq + Nq*(ixco + nxcoDB*(ixnh3 + nxnh3DB*(ixch4 + nxch4DB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(iT)))))); // Working! ! !
                                    //kq[index] = value;
                                    //cout << "Index " << index << " value is " << value << endl;
                                }
                                mixtureStates++;                                
                            }
                        }                    
                    }
                }
            }
        }
        file.close();
        zeroFile.close();
        cout << "Mixture states loaded: " << mixtureStates << ", Overflow states skipped: " << overflowStates << endl;
    }
    cout << "FSK Table (Ammonia Database) successfully loaded into memory." << endl;
    // Go over the kq vector to check for NaN or Inf values
    /*for (size_t i = 0; i < kq.size(); ++i)
    {
        if (std::isnan(kq[i]) || std::isinf(kq[i]))
        {
            //cout << "Warning: kq value is NaN or Inf at index " << i << endl;
            //cout << "Ignoring this value." << endl;
            kq[i] = 0.0f; // Set to zero or handle as needed
        }
        // Fudge number by x% to check for verification (NEEDS TO BE REMOVED)
        //kq[i] = kq[i]*1.10;
    }*/
}

void loadIntoMemFuel(std::vector<float>& kq, std::string fpath, const double Tref)
{
    // FSK Table Parameters
    const int NDB = 32, Nq=32;
    const int nTgDB=28, nxnh3DB = 13, nxch4DB = 13;
    const double xNH3DB[nxnh3DB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xCH4DB[nxch4DB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    int nT=nTgDB;
    // Loading only fuels
    kq.assign(nxch4DB*nxnh3DB*nTgDB*nT*Nq, 0);

    cout << "Loading FSK Table (Fuel DataBase) data into memory." << endl;
    double P = 1.0;

    for (int iTref = 0; iTref < nTgDB; ++iTref) // Tref Loop
    {
        double Tref_loop = 300.0 + iTref * 100.0;
        int TrefInt = static_cast<int>(std::round(Tref_loop));
        //int TrefInt = static_cast<int>(std::round(Tref));
        std::string filepath = fpath + "fsck_Tref_" + std::to_string(TrefInt) + ".bin";
        //cout << "Loading file: " << filepath << endl;
        fstream file;
        file.open(filepath, ios::in | ios::binary);
        if (!file)
        {
            cout << "Error opening file! " << filepath << endl;
            exit(EXIT_FAILURE);
        }
        // Open zero concentration state file
        std::string zeroFilePath = fpath + "zero_concentration_state_" + std::to_string(TrefInt) + ".bin";
        std::ifstream zeroFile(fpath + "zero_concentration_state_" + std::to_string(TrefInt) + ".bin", std::ios::binary);
        if (!zeroFile)
        {
            cout << "Error opening zero concentration state file!" << endl;
            exit(EXIT_FAILURE);
        }
        int mixtureStates = 0;
        int overflowStates = 0;
        for (int iT = 0; iT < nT; ++iT) //fv Loop
        {
            for (int ixch4 = 0; ixch4 < nxch4DB; ++ixch4) //CH4 Loop
            {
                for (int ixnh3 = 0; ixnh3 < nxnh3DB; ++ixnh3) //NH3 Loop
                {
                    double xSum = xNH3DB[ixnh3] + xCH4DB[ixch4];
                    if (xSum > 1.0e0)
                    {
                        overflowStates++;
                        continue;
                    }
                    else if (xSum < 1e-10)
                    {
                        std::vector<float> k_zero(Nq, 0.0f);
                    
                        zeroFile.read(reinterpret_cast<char*>(k_zero.data()), Nq * sizeof(float));
                        for (int iq = 0; iq < Nq; ++iq)
                        {
                            // Only load CH4 and NH3 gases
                            int index = iq + Nq*(ixnh3 + nxnh3DB*(ixch4 + nxch4DB*(iT + nT*(iTref))));
                            kq[index] = k_zero[iq];
                        }                               
                        continue;
                    }

                    for (int iq = 0; iq < NDB; iq++)
                    {
                        float value;
                        file.read(reinterpret_cast<char*>(&value), sizeof(float));
                        // Only load CH4 and NH3 gases
                        int index = iq + Nq*(ixnh3 + nxnh3DB*(ixch4 + nxch4DB*(iT + nT*(iTref))));
                        kq[index] = value;
                    }
                    mixtureStates++;                                
                }                    
            }
        }
        file.close();
        zeroFile.close();
        cout << "Mixture states loaded: " << mixtureStates << ", Overflow states skipped: " << overflowStates << endl;
    }
    cout << "FSK Table (Fuel Database) successfully loaded into memory." << endl;
    // Go over the kq vector to check for NaN or Inf values
    for (size_t i = 0; i < kq.size(); ++i)
    {
        if(kq[i]<0.0f)
        {
            kq[i] = 0.0f; // Set to zero or handle as needed
        }
    }
}

void loadIntoMemProduct(std::vector<float>& kq, std::string fpath, const double Tref)
{
    // FSK Table Parameters
    const int NDB = 32, Nq=32;
    const int nTgDB=28, nxco2DB=13, nxh2oDB=13, nxcoDB=6;
    const double xCO2DB[nxco2DB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xH2ODB[nxh2oDB]={0.0e0,0.010e0,0.020e0,0.030e0,0.040e0,0.050e0,0.10e0,0.15e0,0.20e0,0.25e0,0.50e0,0.75e0,1.0e0};
    const double xCODB[nxcoDB]={0.0e0,0.010e0,0.050e0,0.10e0,0.25e0,0.50e0};
    int nT=nTgDB;
    // Loading only FSK Table v3 gases into memory
    kq.assign(nxco2DB*nxh2oDB*nxcoDB*nTgDB*nT*Nq,0);

    cout << "Loading FSK Table (Products DataBase) data into memory." << endl;
    double P = 1.0;

    for (int iTref = 0; iTref < nTgDB; ++iTref) // Tref Loop
    {
        double Tref_loop = 300.0 + iTref * 100.0;
        int TrefInt = static_cast<int>(std::round(Tref_loop));
        std::string filepath = fpath + "fsck_Tref_" + std::to_string(TrefInt) + ".bin";
        fstream file;
        file.open(filepath, ios::in | ios::binary);
        if (!file)
        {
            cout << "Error opening file! " << filepath << endl;
            exit(EXIT_FAILURE);
        }
        // Open zero concentration state file
        std::string zeroFilePath = fpath + "zero_concentration_state_" + std::to_string(TrefInt) + ".bin";
        std::ifstream zeroFile(fpath + "zero_concentration_state_" + std::to_string(TrefInt) + ".bin", std::ios::binary);
        if (!zeroFile)
        {
            cout << "Error opening zero concentration state file!" << endl;
            exit(EXIT_FAILURE);
        }
        int mixtureStates = 0;
        int overflowStates = 0;
        for (int iT = 0; iT < nT; ++iT) //fv Loop
        {
            for (int ixco2 = 0; ixco2 < nxco2DB; ++ixco2) //CO2 Loop
            {
                for (int ixh2o = 0; ixh2o < nxh2oDB; ++ixh2o) //H2O Loop
                {
                    for (int ixco = 0; ixco < nxcoDB; ++ixco) //CO Loop
                    {
                        double xSum = xCO2DB[ixco2] + xH2ODB[ixh2o] + xCODB[ixco];
                        if (xSum > 1.0e0)
                        {
                            overflowStates++;
                            continue;
                        }
                        else if (xSum < 1e-10)
                        {
                            std::vector<float> k_zero(Nq, 0.0f);
                        
                            zeroFile.read(reinterpret_cast<char*>(k_zero.data()), Nq * sizeof(float));
                            for (int iq = 0; iq < Nq; ++iq)
                            {
                                // Only load CO, CO2, H2O gases
                                int index = iq + Nq*(ixco + nxcoDB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(iT + nT*(iTref)))));
                                kq[index] = k_zero[iq];
                            }                               
                            continue;
                        }

                        for (int iq = 0; iq < NDB; iq++)
                        {
                            float value;
                            file.read(reinterpret_cast<char*>(&value), sizeof(float));

                            // Only load CO, CO2, H2O gases
                            int index = iq + Nq*(ixco + nxcoDB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(iT + nT*(iTref)))));
                            kq[index] = value;
                        }
                        mixtureStates++;                                
                    }
                }
            }
        }
        file.close();
        zeroFile.close();
        cout << "Mixture states loaded: " << mixtureStates << ", Overflow states skipped: " << overflowStates << endl;
    }
    cout << "FSK Table (Product Database) successfully loaded into memory." << endl;
    // Go over the kq vector to check for NaN or Inf values
    for (size_t i = 0; i < kq.size(); ++i)
    {
        if(kq[i]<0.0f)
        {
            kq[i] = 0.0f; // Set to zero or handle as needed
        }
    }
}

/*void get_ka(const GasMixInfo& Mix_Info, const double Tref, const std::vector<float>& kq, const std::vector<double>& g, const std::vector<double>& weight, std::vector<double>& k, std::vector<double>& a)
{
    std::vector<double> kll;
    std::vector<double> gql;

    // Get k-values
    k = get_k(Mix_Info,Tref,kq);

    // Get a-values
    kll = get_k(Mix_Info,Mix_Info.T,kq);
    a = afun(g,k,g,kll,gql,weight);
}*/

void get_ka(const GasMixInfo& Mix_Info, const double Tref, const std::vector<float>& kq, const std::vector<double>& g, const std::vector<double>& weight, const bool isV2, std::vector<double>& k, std::vector<double>& a)
{
    std::vector<double> kll;
    std::vector<double> gql;

    // Get k-values
    if (isV2)
    {
        k = get_k_V2(Mix_Info, Tref, kq);
        // For computing a-values
        kll = get_k_V2(Mix_Info, Mix_Info.T, kq);
        // Print k-values for debugging
        /*cout << "k-values for V2 function at Tref " << Tref << " K:" << endl;
        for (size_t i = 0; i < k.size(); ++i)
        {
            cout << "k[" << i << "] = " << k[i] << endl;
        }*/
    }
    else
    {
        k = get_k(Mix_Info, Tref, kq);
        // For computing a-values
        kll = get_k(Mix_Info, Mix_Info.T, kq);
        // Test out mixing function
        //k = get_k_mix(Mix_Info, Tref, g, kq);
        // For computing a-values
        //kll = get_k_mix(Mix_Info, Mix_Info.T, g, kq);
    }
    
    // Get a-values
    a = afun(g,k,g,kll,gql,weight);
}

void get_ka_mix(const GasMixInfo& Mix_Info, const double Tref, const std::vector<float>& kq, const std::vector<float>& kq2, const std::vector<double>& g, const std::vector<double>& weight, std::vector<double>& k, std::vector<double>& a, const int mix_scheme)
{
    std::vector<double> kll;
    std::vector<double> gql;

    // Get k-values
    k = get_k_mix(Mix_Info, Tref, g, kq, kq2, weight, mix_scheme);
    // For computing a-values
    kll = get_k_mix(Mix_Info, Mix_Info.T, g, kq, kq2, weight, mix_scheme);
    
    // Get a-values
    a = afun(g,k,g,kll,gql,weight);
}

void quadgen2(const bool cheb2, std::vector<double>& g, std::vector<double>& weight, const int points, const double alpha)
{
    int n = 2*points;
    std::vector<double> gg(n,0.0);
    std::vector<double> ww(n,0.0);
    const double pi = 3.1415926535897932384626;
    double theta,sum;

    // gaucheb2 snippet
    for (int k = 1; k <= (n+1)/2; ++k)
    {
        theta = double(k)*pi/(n+1.0);
        gg[k-1] = std::cos(theta);
        sum = 0.0;
        for (int m = 1; m <= (n+1)/2; ++m)
        {
            sum = sum + std::sin(((2*m-1)*theta))/(2.0*m-1.0);
        }
        ww[k-1] = 4.0*std::sin(theta)*sum/(n+1.0);
        gg[n+1-k-1] = -gg[k-1];
        ww[n+1-k-1] = ww[k-1];
    }
    // gaucheb2 END

    for (int k = 1; k <= points; ++k)
    {
        g[k-1] = -gg[points+k-1];
        weight[k-1] = ww[points+k-1];
    }
    sum = 0.0;
    if (cheb2 == false && alpha < 0.0)
    {
        for (int i = 1; i <= points; ++i)
        {
            weight[i-1] = 1.5*std::sqrt(1.0-g[i-1])*weight[i-1];
            g[i-1] = 1.0-std::pow((1.0-g[i-1]),1.5);
            sum = sum + weight[i-1];
        }
    }
    else if (cheb2 == false && alpha > 0.0)
    {
        for (int i = 1; i <= points; ++i)
        {
            weight[i-1] = alpha*weight[i-1]*std::pow((1.0-g[i-1]),alpha-1.0);
            g[i-1] = 1.0-std::pow((1.0-g[i-1]),alpha);
            sum = sum + weight[i-1];
        }
    }
    for (int i = 1; i <= points; ++i)
    {
        weight[i-1] = weight[i-1]/sum;
    }
}

std::vector<double> get_k(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq)
{
    const int Nq=32; //, NDB = 32, Ndm=7;
    const int nTgDB=28, nxco2DB=13, nxh2oDB=13, nxcoDB=6; //, nfvDB=2;
    int nT = nTgDB;

    std::vector<std::vector<int>> grid;
    std::vector<int> lp;
    std::vector<double> coef;
    std::vector<double> k(Nq,0.0);  // Output vector with k values

    coef_setup(Mix_Info,Tref,grid,coef);
    lp = detmLP(coef);

    for (int iP = 1; iP<=lp[0]; ++iP)
    {
        // Pressure loop not functional for Table V3
        // Only 1 atm pressure available but the loop is included
        // as per source.
        double wP = (iP==1) ? 1.0-coef[0]:coef[0];
        for (int ifv = 1; ifv <= lp[1]; ++ifv)
        {
            double wfv = (ifv==1) ? 1.0-coef[1]:coef[1];
            for (int ixco2 = 1; ixco2 <= lp[2]; ++ixco2)
            {
                double wxco2 = (ixco2==1) ? 1.0-coef[2]:coef[2];
                for (int ixh2o = 1; ixh2o <= lp[3]; ++ixh2o)
                {
                    double wxh2o = (ixh2o==1) ? 1.0-coef[3]:coef[3];
                    for (int ixco = 1; ixco <= lp[4]; ++ixco)
                    {
                        double wxco = (ixco==1) ? 1.0-coef[4]:coef[4];
                        for (int iTg = 1; iTg <= lp[5]; ++iTg)
                        {
                            double wTg = (iTg==1) ? 1.0-coef[5]:coef[5];
                            for (int iT = 1; iT <= lp[6]; ++iT)
                            {
                                double wT = (iT==1) ? 1.0-coef[6]:coef[6];

                                double we = wP*wfv*wxco2*wxh2o*wxco*wTg*wT;

                                for (int iq = 0; iq < Nq; ++iq)
                                {
                                    int index = iq + Nq*(grid[6][iT-1] + nT*(grid[5][iTg-1] + nTgDB*(grid[4][ixco-1] + nxcoDB*(grid[3][ixh2o-1] + nxh2oDB*(grid[2][ixco2-1] + nxco2DB*(grid[1][ifv-1]))))));
                                    //int index = iq + Nq*(iT + nT*(iTg + nTgDB*(ixco + nxcoDB*(ixh2o + nxh2oDB*(ixco2 + nxco2DB*(ifv))))));
                                    k[iq] = k[iq] + we*kq[index];
                                } // Quadrature point loop
                            } // iT loop
                        } // iTg loop
                    } // ixco loop
                } // ixh2o loop
            } // ixco2 loop
        } // ifv loop
    } // iP loop

    return k;
}

std::vector<double> get_k_V2(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq)
{
    const int Nq=32; //, NDB = 32, Ndm=7;
    const int nTgDB=28, nxco2DB=13, nxh2oDB=13, nxcoDB=6, nxnh3DB = 13, nxch4DB = 13; //, nfvDB=2;
    int nT = nTgDB;

    std::vector<std::vector<int>> grid;
    std::vector<int> lp;
    std::vector<double> coef;
    std::vector<double> k(Nq,0.0);  // Output vector with k values

    coef_setupV2(Mix_Info,Tref,grid,coef);
    lp = detmLP(coef);

    for (int iTg = 1; iTg <= lp[8]; ++iTg)
    {
        double wTg = (iTg==1) ? 1.0-coef[8]:coef[8];
        for (int iT = 1; iT <= lp[7]; ++iT)
        {
            double wT = (iT==1) ? 1.0-coef[7]:coef[7];
            for (int ixco2 = 1; ixco2 <= lp[2]; ++ixco2)
            {
                double wxco2 = (ixco2==1) ? 1.0-coef[2]:coef[2];
                for (int ixh2o = 1; ixh2o <= lp[3]; ++ixh2o)
                {
                    double wxh2o = (ixh2o==1) ? 1.0-coef[3]:coef[3];
                    for (int ixch4 = 1; ixch4 <= lp[4]; ++ixch4)
                    {
                        double wxch4 = (ixch4==1) ? 1.0-coef[4]:coef[4];
                        for (int ixnh3 = 1; ixnh3 <= lp[5]; ++ixnh3)
                        {
                            double wxnh3 = (ixnh3==1) ? 1.0-coef[5]:coef[5];
                            for (int ixco = 1; ixco <= lp[6]; ++ixco)
                            {
                                double wxco = (ixco==1) ? 1.0-coef[6]:coef[6];
                                
                                //double we = wxco2*wxh2o*wxch4*wxnh3*wxco*wT*wTg;
                                //double we = wxco2*wxh2o*wxco*wTg*wT;
                                double we = wxnh3*wxch4*wTg*wT;

                                for (int iq = 0; iq < Nq; ++iq)
                                {
                                    // Only load CO, CO2, H2O gases loaded
                                    //int index = iq + Nq*(grid[6][ixco-1] + nxcoDB*(grid[3][ixh2o-1] + nxh2oDB*(grid[2][ixco2-1] + nxco2DB*(grid[7][iT-1] + nT*(grid[8][iTg-1])))));
                                    // Only load CH4 and NH3 gases loaded
                                    int index = iq + Nq*(grid[5][ixnh3-1] + nxnh3DB*(grid[4][ixch4-1] + nxch4DB*(grid[7][iT-1] + nT*(grid[8][iTg-1]))));
                                    
                                    // All gases loaded
                                    //int index = iq + Nq*(grid[6][ixco-1] + nxcoDB*(grid[5][ixnh3-1] + nxnh3DB*(grid[4][ixch4-1] + nxch4DB*(grid[3][ixh2o-1] + nxh2oDB*(grid[2][ixco2-1] + nxco2DB*(grid[7][iT-1]))))));
                                    k[iq] = k[iq] + we*kq[index];
                                    //cout << "Index = " << index << ", kq[index] = " << kq[index] << endl;
                                } // Quadrature point loop                            
                            } // ixco loop
                        } // ixnh3 loop
                    } // ixch4 loop
                } // ixh2o loop
            } // ixco2 loop
        } // iT loop
    } // iTg loop
    return k;
}

std::vector<double> get_k_fuel(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq)
{
    const int Nq=32; //, NDB = 32, Ndm=7;
    const int nTgDB=28, nxnh3DB = 13, nxch4DB = 13; //, nfvDB=2;
    int nT = nTgDB;

    std::vector<std::vector<int>> grid;
    std::vector<int> lp;
    std::vector<double> coef;
    std::vector<double> k(Nq,0.0);  // Output vector with k values

    coef_setupV2(Mix_Info,Tref,grid,coef);
    lp = detmLP(coef);

    for (int iTg = 1; iTg <= lp[8]; ++iTg)
    {
        double wTg = (iTg==1) ? 1.0-coef[8]:coef[8];
        for (int iT = 1; iT <= lp[7]; ++iT)
        {
            double wT = (iT==1) ? 1.0-coef[7]:coef[7];
            for (int ixch4 = 1; ixch4 <= lp[4]; ++ixch4)
            {
                double wxch4 = (ixch4==1) ? 1.0-coef[4]:coef[4];
                for (int ixnh3 = 1; ixnh3 <= lp[5]; ++ixnh3)
                {
                    double wxnh3 = (ixnh3==1) ? 1.0-coef[5]:coef[5];
                    // Calculate totatl weight
                    double we = wxnh3*wxch4*wTg*wT;

                    for (int iq = 0; iq < Nq; ++iq)
                    {
                        // Only load CH4 and NH3 gases loaded
                        int index = iq + Nq*(grid[5][ixnh3-1] + nxnh3DB*(grid[4][ixch4-1] + nxch4DB*(grid[7][iT-1] + nT*(grid[8][iTg-1]))));
                        k[iq] = k[iq] + we*kq[index];
                    } // Quadrature point loop                            
                } // ixnh3 loop
            } // ixch4 loop
        } // iT loop
    } // iTg loop
    return k;
}

std::vector<double> get_k_product(const GasMixInfo& Mix_Info, double Tref, const std::vector<float>& kq)
{
    const int Nq=32; //, NDB = 32, Ndm=7;
    const int nTgDB=28, nxco2DB=13, nxh2oDB=13, nxcoDB=6;
    int nT = nTgDB;

    std::vector<std::vector<int>> grid;
    std::vector<int> lp;
    std::vector<double> coef;
    std::vector<double> k(Nq,0.0);  // Output vector with k values

    coef_setupV2(Mix_Info,Tref,grid,coef);
    lp = detmLP(coef);

    for (int iTg = 1; iTg <= lp[8]; ++iTg)
    {
        double wTg = (iTg==1) ? 1.0-coef[8]:coef[8];
        for (int iT = 1; iT <= lp[7]; ++iT)
        {
            double wT = (iT==1) ? 1.0-coef[7]:coef[7];
            for (int ixco2 = 1; ixco2 <= lp[2]; ++ixco2)
            {
                double wxco2 = (ixco2==1) ? 1.0-coef[2]:coef[2];
                for (int ixh2o = 1; ixh2o <= lp[3]; ++ixh2o)
                {
                    double wxh2o = (ixh2o==1) ? 1.0-coef[3]:coef[3];
                    for (int ixco = 1; ixco <= lp[6]; ++ixco)
                    {
                        double wxco = (ixco==1) ? 1.0-coef[6]:coef[6];
                        // Calculate total weight
                        double we = wxco2*wxh2o*wxco*wTg*wT;

                        for (int iq = 0; iq < Nq; ++iq)
                        {
                            // Only load CO, CO2, H2O gases loaded
                            int index = iq + Nq*(grid[6][ixco-1] + nxcoDB*(grid[3][ixh2o-1] + nxh2oDB*(grid[2][ixco2-1] + nxco2DB*(grid[7][iT-1] + nT*(grid[8][iTg-1])))));
                            k[iq] = k[iq] + we*kq[index];
                        } // Quadrature point loop                            
                    } // ixco loop
                } // ixh2o loop
            } // ixco2 loop
        } // iT loop
    } // iTg loop
    return k;
}

std::vector<double> afun(const std::vector<double> grr,const std::vector<double> krr,const std::vector<double> grc,const std::vector<double> krc,const std::vector<double> g, const std::vector<double> w)
{
    int nq = w.size();
    std::vector<double> a(nq,0.0);

    std::vector<double> gb(nq+1,0.0);
    std::vector<double> gq(nq,0.0);
    std::vector<double> kq(nq,0.0);

    for (int iq = 1; iq <= nq; ++iq)
    {
        gb[iq] = gb[iq-1]+w[iq-1];
    }

    kq = linearInterpMono(grr,krr,std::vector<double>(gb.begin()+1,gb.end()));
    gq = linearInterpMono(krc,grc,kq);
    gq.insert(gq.begin(),0.0);
    if (static_cast<int>(gq.size()) != nq+1)
    {
        cout << "Check code for some error!" << endl;
    }
    for (int i = 1; i <= nq; ++i)
    {
        a[i-1] = (gq[i]-gq[i-1])/w[i-1];
    }

    return a;
}

std::vector<double> linearInterpMono(const std::vector<double> xx, const std::vector<double> yy, const std::vector<double> xi)
{
    int nxy = static_cast<int>(xx.size());
    int ni = static_cast<int>(xi.size());
    std::vector<double> yi(ni,0.0); // Return vector

    int n = 1;
    for (int i = 0; i < ni; ++i)
    {
        while (xi[i] >= xx[n-1])
        {
            n = n+1;
            if(n>nxy)
            {
                std::fill(yi.begin()+i,yi.end(),yy[nxy-1]);
                return yi;
            }
        }
        if (n==1)
        {
            yi[i] = yy[0];
        }
        else
        {
            yi[i] = yy[n-1]+(yy[n-2]-yy[n-1])*(xi[i]-xx[n-1])/(xx[n-2]-xx[n-1]);
        }
    }
    return yi;
}

void coef_setup(const GasMixInfo& Mix_Info, const double Tref, std::vector<std::vector<int>>& grid, std::vector<double>& coef)
{
    const int Ndm = 7;
    const int nP =1;
    const int nTgDB=28, nxco2=13, nxh2o=13, nxco=6, nfv=2;
    int nT=nTgDB; int nTg=nT;
    const double dTg = 100.0;
    const double dT = 100.0;
    //const double xCO2DB[13]={0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    //const double xH2ODB[13]={0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    //const double xCODB[6]={0.e0,0.01e0,0.05e0,0.1e0,0.25e0,0.5e0};
    //const double fvDB[2]={0.e0,1e-5};
    std::vector<double> xCO2 = {0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    std::vector<double> xH2O = {0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    std::vector<double> xCO = {0.e0,0.01e0,0.05e0,0.1e0,0.25e0,0.5e0};
    std::vector<double> fv = {0.e0,1e-5};
    std::vector<double> P(nP,0.0);
    std::vector<double> T(nT,0.0);
    std::vector<double> Tg;

    P[1]=1.0;
    for (int i = 0; i<nT; ++i)
    {
        T[i]=300.0+i*dT;
    }
    Tg = T;

    grid.assign(Ndm, std::vector<int>(2,0));
    coef.assign(Ndm, 0.0);

    // Pressure Index
    if (nP == 1)
    {
        grid[0][0]=0;
        coef[0]=0.0;
    }
    else
    {
        grid[0][0]=detmID(P,Mix_Info.P);
        grid[0][1] = grid[0][0] + 1;
        coef[0] = (Mix_Info.P - P[grid[0][0]]) / (P[grid[0][1]] - P[grid[0][0]]);
    }

    // fv Index
    if (nfv == 1)
    {
        grid[1][0]=0;
        coef[1]=0.0;
    }
    else
    {
        grid[1][0]=detmID(fv,Mix_Info.fv);
        grid[1][1] = grid[1][0] + 1;
        coef[1] = (Mix_Info.fv - fv[grid[1][0]]) / (fv[grid[1][1]] - fv[grid[1][0]]);
    }

    // CO2 Index
    if (nxco2 == 1)
    {
        grid[2][0]=0;
        coef[2]=0.0;
    }
    else
    {
        grid[2][0]=detmID(xCO2,Mix_Info.xCO2);
        grid[2][1] = grid[2][0] + 1;
        coef[2] = (Mix_Info.xCO2 - xCO2[grid[2][0]]) / (xCO2[grid[2][1]] - xCO2[grid[2][0]]);
    }

    // H2O Index
    if (nxh2o == 1)
    {
        grid[3][0]=0;
        coef[3]=0.0;
    }
    else
    {
        grid[3][0]=detmID(xH2O,Mix_Info.xH2O);
        grid[3][1] = grid[3][0] + 1;
        coef[3] = (Mix_Info.xH2O - xH2O[grid[3][0]]) / (xH2O[grid[3][1]] - xH2O[grid[3][0]]);
    }

    // CO Index
    if (nxco == 1)
    {
        grid[4][0]=0;
        coef[4]=0.0;
    }
    else
    {
        grid[4][0]=detmID(xCO,Mix_Info.xCO);
        grid[4][1] = grid[4][0] + 1;
        coef[4] = (Mix_Info.xCO - xCO[grid[4][0]]) / (xCO[grid[4][1]] - xCO[grid[4][0]]);
    }

    // Local Temperature index
    grid[5][0] = static_cast<int>((Mix_Info.T - Tg[0]) / dTg);
    grid[5][0] = max(grid[5][0],0);
    grid[5][0] = min(grid[5][0],nTg - 1 - 1);
    grid[5][1] = grid[5][0] + 1;
    coef[5] = (Mix_Info.T - Tg[grid[5][0]]) / dTg;

    // Reference Temperature index
    grid[6][0] = static_cast<int>((Tref - T[0]) / dT);
    grid[6][0] = max(grid[6][0],0);
    grid[6][0] = min(grid[6][0],nT - 1 - 1);
    grid[6][1] = grid[6][0] + 1;
    coef[6] = (Tref - T[grid[6][0]]) / dT;
}

void coef_setupV2(const GasMixInfo& Mix_Info, const double Tref, std::vector<std::vector<int>>& grid, std::vector<double>& coef)
{
    const int Ndm = 9;
    const int nP =1;
    const int nTgDB=28, nxco2=13, nxh2o=13, nxco=6, nxnh3=13, nxch4=13; //, nfv=2;
    int nT=nTgDB; int nTg=nT;
    const double dTg = 100.0;
    const double dT = 100.0;

    std::vector<double> xCO2 = {0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    std::vector<double> xH2O = {0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    std::vector<double> xNH3 = {0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    std::vector<double> xCH4 = {0.e0,0.01e0,0.02e0,0.03e0,0.04e0,0.05e0,0.1e0,0.15e0,0.2e0,0.25e0,0.5e0,0.75e0,1.e0};
    std::vector<double> xCO = {0.e0,0.01e0,0.05e0,0.1e0,0.25e0,0.5e0};
    //std::vector<double> fv = {0.e0,1e-5};
    std::vector<double> P(nP,0.0);
    std::vector<double> T(nT,0.0);
    std::vector<double> Tg;
    P[1]=1.0;
    for (int i = 0; i<nT; ++i)
    {
        T[i]=300.0+i*dT;
    }
    Tg = T;

    grid.assign(Ndm, std::vector<int>(2,0));
    coef.assign(Ndm, 0.0);

    // Pressure Index
    // Fill in the pressure index when variety of pressures is available
    // grid[0] variables must be used

    // Soot volume fraction Index
    // Fill in when variety of soot volume fractions is available (nfv > 1) 
    // grid[1] variables must be used

    // CO2 Index
    if (nxco2 == 1)
    {
        grid[2][0]=0;
        coef[2]=0.0;
    }
    else
    {
        grid[2][0]=detmID(xCO2,Mix_Info.xCO2);
        grid[2][1] = grid[2][0] + 1;
        coef[2] = (Mix_Info.xCO2 - xCO2[grid[2][0]]) / (xCO2[grid[2][1]] - xCO2[grid[2][0]]);
    }

    // H2O Index
    if (nxh2o == 1)
    {
        grid[3][0]=0;
        coef[3]=0.0;
    }
    else
    {
        grid[3][0]=detmID(xH2O,Mix_Info.xH2O);
        grid[3][1] = grid[3][0] + 1;
        coef[3] = (Mix_Info.xH2O - xH2O[grid[3][0]]) / (xH2O[grid[3][1]] - xH2O[grid[3][0]]);
    }

    // CH4 Index
    if (nxch4 == 1)
    {
        grid[4][0]=0;
        coef[4]=0.0;
    }
    else
    {
        grid[4][0]=detmID(xCH4,Mix_Info.xCH4);
        grid[4][1] = grid[4][0] + 1;
        coef[4] = (Mix_Info.xCH4 - xCH4[grid[4][0]]) / (xCH4[grid[4][1]] - xCH4[grid[4][0]]);
    }

    // NH3 Index
    if (nxnh3 == 1)
    {
        grid[5][0]=0;
        coef[5]=0.0;
    }
    else
    {
        grid[5][0]=detmID(xNH3,Mix_Info.xNH3);
        grid[5][1] = grid[5][0] + 1;
        coef[5] = (Mix_Info.xNH3 - xNH3[grid[5][0]]) / (xNH3[grid[5][1]] - xNH3[grid[5][0]]);
    }

    // CO Index
    if (nxco == 1)
    {
        grid[6][0]=0;
        coef[6]=0.0;
    }
    else
    {
        grid[6][0]=detmID(xCO,Mix_Info.xCO);
        grid[6][1] = grid[6][0] + 1;
        coef[6] = (Mix_Info.xCO - xCO[grid[6][0]]) / (xCO[grid[6][1]] - xCO[grid[6][0]]);
    }

    // Local Temperature index
    grid[7][0] = static_cast<int>((Mix_Info.T - Tg[0]) / dTg);
    grid[7][0] = max(grid[7][0],0);
    grid[7][0] = min(grid[7][0],nTg - 1 - 1);
    grid[7][1] = grid[7][0] + 1;
    coef[7] = (Mix_Info.T - Tg[grid[7][0]]) / dTg;

    // Reference Temperature Index (NEEDED FOR STRETCHING FACTORS!)
    grid[8][0] = static_cast<int>((Tref - T[0]) / dT);
    grid[8][0] = max(grid[8][0],0);
    grid[8][0] = min(grid[8][0],nT - 1 - 1);
    grid[8][1] = grid[8][0] + 1;
    coef[8] = (Tref - T[grid[8][0]]) / dT;
}

std::vector<int> detmLP(std::vector<double> arr)
{
    //size_t n = sizeof(arr)/sizeof(double);
    //int* val = new(nothrow) int[n];
    int n = static_cast<int>(arr.size());
    std::vector<int> val(n,0);
    double cmin = 1e-3;
    for (int i = 0; i < n; i++)
    {
        if (arr[i] <= cmin)
        {
            val[i]=1;
        }
        else
        {
            val[i]=2;
        }
    }
    return val;
}

int detmID(std::vector<double> arr, double val)
{
    // CHECK CODE AFTER ALTERING TO USE VECTORS
    int id;
    int lef, rig, mid;

    if (val >= arr[arr.size()-1])
    {
        id = static_cast<int>(arr.size()) - 1 - 1;
    }
    else if (val <= arr[0])
    {
        id = 0;
    }
    else
    {
        lef = 0;
        rig = static_cast<int>(arr.size()) - 1;
        mid = (lef+rig)/2;
        while (rig > lef +1)
        {
            if (val < arr[mid])
            {
                rig = mid;
            }
            else
            {
                lef = mid;
            }
            mid = (lef+rig)/2;
        }
        id = lef;
    }
    return id;
}

// Functions for mixing schemes
std::vector<double> get_k_mix(const GasMixInfo& Mix_Info, double Tref, const std::vector<double>& g, const std::vector<float>& kq, const std::vector<float>& kq2, const std::vector<double>& weight, const int mix_scheme)
{
    std::vector<double> k1, k2;
    std::vector<double> kmix, gmix, kmix_quad;

    // Get k-values for individual gases
    // Products of combustion table
    GasMixInfo Mix_Info1 = Mix_Info;
    Mix_Info1.xNH3 = 0.0; Mix_Info1.xCH4 = 0.0; 
    //k1 = get_k(Mix_Info1, Tref, kq);
    k1 = get_k_product(Mix_Info1, Tref, kq);
    // Fuel of combustion table
    GasMixInfo Mix_Info2 = Mix_Info;
    Mix_Info2.xCO2 = 0.0; Mix_Info2.xH2O = 0.0; Mix_Info2.xCO = 0.0;
    //k2 = get_k_V2(Mix_Info2, Tref, kq2);
    k2 = get_k_fuel(Mix_Info2, Tref, kq2);
    
    // Mix k-values
    double kmin = k1.front() + k2.front();
    double kmax = k1.back() + k2.back();
    kmix = kPowerLaw(kmin, kmax, 128, 0.1e0);

    // Select the mixing scheme
    switch (mix_scheme)
    {
    case 1: // Multiplication
        gmix = multiplication(k1, g, k2, g, kmix);
        break;
    case 2: // Superposition
        gmix = superposition(k1, g, k2, g, kmix);
        break;
    case 3: // Modest Riazzi
        gmix = modestRiazzi(k1, g, k2, g, kmix, weight);
        break;
    default:
        std::cerr << "Invalid mixing scheme selected. Using multiplication by default." << std::endl;
        gmix = multiplication(k1, g, k2, g, kmix);
        break;
    }
    //gmix = multiplication(k1, g, k2, g, kmix);
    //gmix = superposition(k1, g, k2, g, kmix);
    //gmix = modestRiazzi(k1, g, k2, g, kmix, weight);

    // Write out the mixed k-values and g-values
    /*cout << "Mixed k-values:" << endl;
    for (size_t i = 0; i < kmix.size(); ++i) 
    {
        cout << "k[" << i << "] = " << kmix[i] << ", g[" << i << "] = " << gmix[i] << endl;
    }*/

    // Interpolate to the quadrature points
    kmix_quad = linearInterpMono(gmix, kmix, g);

    return kmix_quad;
}

std::vector<double> multiplication(
    const std::vector<double>& k1,
    const std::vector<double>& g1,
    const std::vector<double>& k2,
    const std::vector<double>& g2,
    const std::vector<double>& k)
{
    int nk1 = static_cast<int>(k1.size());
    int nk2 = static_cast<int>(k2.size());
    int nop = static_cast<int>(k.size());

    std::vector<double> g(nop), gq1(nop), gq2(nop);

    //linearInterpMono(nk1, k1, g1, nop, k, gq1);
    //linearInterpMono(nk2, k2, g2, nop, k, gq2);

    gq1 = linearInterpMono(k1, g1, k);
    gq2 = linearInterpMono(k2, g2, k);

    for (int i = 0; i < nop; ++i) {
        g[i] = gq1[i] * gq2[i];
    }

    double kmin = k1.front() + k2.front();
    double kmax = k1.back() + k2.back();
    bound(k, g, kmin, kmax);

    return g;
}

std::vector<double> superposition(
    const std::vector<double>& k1,
    const std::vector<double>& g1,
    const std::vector<double>& k2,
    const std::vector<double>& g2,
    const std::vector<double>& k)
{
    int nk1 = static_cast<int>(k1.size());
    int nk2 = static_cast<int>(k2.size());
    int nop = static_cast<int>(k.size());

    std::vector<double> g(nop), gq1(nop), gq2(nop);

    gq1 = linearInterpMono(k1, g1, k);
    gq2 = linearInterpMono(k2, g2, k);

    for (int i = 0; i < nop; ++i) {
        g[i] = gq1[i] + gq2[i] - 1.0e0;
    }

    double kmin = k1.front() + k2.front();
    double kmax = k1.back() + k2.back();
    bound(k, g, kmin, kmax);

    return g;
}

std::vector<double> modestRiazzi(
    const std::vector<double>& k1,
    const std::vector<double>& g1,
    const std::vector<double>& k2,
    const std::vector<double>& g2,
    const std::vector<double>& k,
    const std::vector<double>& wq)
{
    int nk1 = static_cast<int>(k1.size());
    int nk2 = static_cast<int>(k2.size());
    int nop = static_cast<int>(k.size());

    std::vector<double> g(nop, 0.0);

    double kmin = k1.front() + k2.front();
    double kmax = k1.back() + k2.back();

    for (int i = 0; i < static_cast<int>(k.size()); ++i) 
    {
        if (k[i] < kmin) 
        {
            g[i] = 0.0;
            continue; // k is too small, next loop
        }
        if (k[i] > kmax) 
        {
            break; // k is too large, no need to loop other k's
        }

        // Compute kq2 as k[i] - kq1 (element-wise)
        std::vector<double> kq2(k1.size());
        for (size_t j = 0; j < k1.size(); ++j) 
        {
            kq2[j] = k[i] - k1[j];
        }

        // Reverse kq2 for descending order
        std::vector<double> kq2_rev(kq2.rbegin(), kq2.rend());

        // Interpolate: gq2_rev = linearInterpMono(k2, g2, kq2_rev)
        std::vector<double> gq2_rev = linearInterpMono(k2, g2, kq2_rev);

        // Reverse gq2_rev to get gq2 in original order
        std::vector<double> gq2(gq2_rev.rbegin(), gq2_rev.rend());

        // Lower bound zero
        for (double& val : gq2) 
        {
            val = std::max(val, 0.0);
        }

        // Weighted sum
        double sum = 0.0;
        for (size_t j = 0; j < gq2.size(); ++j) 
        {
            sum += gq2[j] * wq[j];
        }
        g[i] = sum;
    }

    return g;
}

void bound(const std::vector<double>& k, std::vector<double>& g, double kmin, double kmax)
{
    const double ZERO = 0.0e0, ONE = 1.0e0;
    int nop = static_cast<int>(k.size());

    // Clamp g between 0 and 1
    for (int i = 0; i < nop; ++i) {
        g[i] = std::max(g[i], ZERO);
        g[i] = std::min(g[i], ONE);
    }

    // Apply kmin/kmax logic
    for (int i = 0; i < nop; ++i) {
        if (k[i] < kmin) {
            g[i] = ZERO;
        } else if (k[i] > kmax) {
            g[i] = ONE;
        }
    }
}

std::vector<double> kPowerLaw(double kmin, double kmax, int n, double pwr)
{
    std::vector<double> k(n, 0.0);
    double pwrk_min = std::pow(kmin, pwr);
    double pwrk_max = std::pow(kmax, pwr);
    double pwrk_step = (pwrk_max - pwrk_min) / static_cast<double>(n - 1);

    for (int i = 0; i < n; ++i)
    {
        double val = pwrk_min + static_cast<double>(i) * pwrk_step;
        k[i] = std::pow(val, 1.0 / pwr);
    }
    return k;
}