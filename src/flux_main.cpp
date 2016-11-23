#include <wipp.h>

using namespace std;
using namespace Eigen;

#pragma omp 
int main(int argc, char *argv[]) 
{
    FILE* inPtr;
    string inp_dir;
    string out_dir;
    string flux_filename;
    ostringstream alpha_N_name, alpha_S_name;

    vector <double> freqs;
    double out_lat, out_lon;

    double dA_N[NUM_E][NUM_TIMES] = {0};
    double dA_S[NUM_E][NUM_TIMES] = {0};


    // Default parameters:
    inp_dir = "/shared/users/asousa/WIPP/3dWIPP/outputs";
    out_dir = "/shared/users/asousa/WIPP/3dWIPP/outputs";
    flux_filename = "/shared/users/asousa/WIPP/3dWIPP/data/EQFLUXMA.dat";
    out_lat = 45;
    out_lon = 0;



    // Parse input arguments:
    static struct option long_options[] =
    {
        {"inp_dir",     required_argument,    0, 'a'},
        {"flux_file",   required_argument,    0, 'b'},
        {"lat",         required_argument,    0, 'c'},
        {"lon",         required_argument,    0, 'd'},
        {0, 0, 0, 0}
    };

    int opt = 0;
    int opt_index = 0;

    while (opt != -1) {
        opt = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:", long_options, &opt_index);
        // cout << "opt is " << opt << "\n";
        switch(opt) {
            case 0:
            if (long_options[opt_index].flag != 0)      break;
            case 'a':   // inp_dir
                inp_dir = (string) optarg;              break;
            case 'b':
                flux_filename = (string) optarg;        break;
            case 'c':
                out_lat = strtod(optarg, NULL);         break;
            case 'd':
                out_lon = strtod(optarg, NULL);         break;
            case 'f':
                freqs.push_back(strtod(optarg,NULL));   break;
            case '?':
                 printf("\nUnknown option: %s\n",opt);
            break;
        }
    }


    cout << "---- 3D WIPP ----\n";
    cout << "   Flux Calculation\n";
    cout << "-----------------\n\n";

    cout << "---- Input Parameters ----\n";
    cout << "Input dir:\t\t" << inp_dir <<"\n";
    cout << "Flux file:\t\t" << flux_filename <<"\n";
    cout << "Latitude:\t\t"  << out_lat << "\n";
    cout << "longitude:\t\t" << out_lon << "\n";
    cout << "Frequencies:\t\t";
    print_vector(freqs);

    // Output file names:
    alpha_N_name << out_dir << "/alpha_" << out_lat << "_" << out_lon << "N.dat";
    alpha_S_name << out_dir << "/alpha_" << out_lat << "_" << out_lon << "S.dat";

    // Load all p-files and sum them:
    
    for (int k=0; k<1; ++k) {
        string NS = (k==0 ? "N" : "S");
        for ( int f_ind=0; f_ind < freqs.size(); ++f_ind) {
            ostringstream pname;
            pname << inp_dir << "/p" << NS << "_" << out_lat << "_" << out_lon << "_" << freqs[f_ind] << ".dat";

            double ptmp[NUM_E][NUM_TIMES];

            read_p_array(ptmp, pname.str());
            
            for (int row=0; row < NUM_E; ++row) {
                for (int col=0; col < NUM_TIMES; ++col){        
                    if (k==0) {
                        dA_N[row][col]+= ptmp[row][col];
                    } else {
                        dA_S[row][col]+= ptmp[row][col];
                    }
                }
            }
        }
    }
}