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
    ostringstream summed_N_name, summed_S_name;
    int alpha_dist, flux_dist;
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
    alpha_dist = 0; // Which pitch-angle distribution to use: 0=ramp, 1=square 
    flux_dist = 0;  // Which flux distribution to use: 0=ae8 file, 1=Bell 2002, 2=flat.



    // Parse input arguments:
    static struct option long_options[] =
    {
        {"inp_dir",     required_argument,    0, 'a'},
        {"out_dir",     required_argument,    0, 'b'},
        {"flux_file",   required_argument,    0, 'c'},
        {"lat",         required_argument,    0, 'd'},
        {"lon",         required_argument,    0, 'e'},
        {"f",           required_argument,    0, 'f'},
        {"alpha_dist",  required_argument,    0, 'g'},
        {"flux_dist",   required_argument,    0, 'h'},
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
                out_dir = (string) optarg;              break;
            case 'c':
                flux_filename = (string) optarg;        break;
            case 'd':
                out_lat = strtod(optarg, NULL);         break;
            case 'e':
                out_lon = strtod(optarg, NULL);         break;
            case 'f':
                freqs.push_back(strtod(optarg,NULL));   break;
            case 'g':
                alpha_dist = strtod(optarg, NULL);      break;
            case 'h':
                flux_dist = strtod(optarg, NULL);       break;
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
    cout << "Frequencies:\t\t"; print_vector(freqs);
    cout << "Flux distribution:\t\t" << flux_dist << "\n";
    cout << "Alpha distribution:\t\t" << alpha_dist << "\n";

    // Output file names:
    alpha_N_name << out_dir << "/alpha_" << out_lat << "_" << out_lon << "N.dat";
    alpha_S_name << out_dir << "/alpha_" << out_lat << "_" << out_lon << "S.dat";


    // // Load flux file:
    double J[100][100];
    readJ(J, flux_filename);

    // Print out the J array for debuggins
    // for (int x=0; x<100; ++x) {
    //     cout << x << ": ";
    //     for (int y=0; y<100; ++y) {
    //         cout << J[x][y] << " ";
    //     }
    //     cout << endl;
    // }

    // Load all p-files and sum them:
    for (int k=0; k<=1; ++k) {

        string NS = (k==0 ? "N" : "S");
        cout << NS << endl;
        for ( int f_ind=0; f_ind < freqs.size(); ++f_ind) {
            ostringstream pname, pname_zipped, zip_cmd;

            pname << inp_dir << "/p" << NS << "_" << out_lat << "_" << out_lon << "_" << freqs[f_ind] << ".dat";
            zip_cmd << "gunzip " << pname.str() << ".gz";
            double ptmp[NUM_E][NUM_TIMES];

            // Unzip
            exec(zip_cmd.str().c_str());
            // Read
            read_p_array(ptmp, pname.str());

            // Sum
            for (int row=0; row < NUM_E; ++row) {
                for (int col=0; col < NUM_TIMES; ++col){
                    if (k==0) {
                        dA_N[row][col]+= ptmp[row][col];
                    } else {
                        dA_S[row][col]+= ptmp[row][col];
                    }
                }
            }

            // Re-zip
            zip_cmd.str("");
            zip_cmd << "gzip " << pname.str();
            exec(zip_cmd.str().c_str());
        }
    }

    // Save the summed arrays:
    summed_N_name << out_dir << "/Psum_" << out_lat << "_" << out_lon << "_N.dat";
    summed_S_name << out_dir << "/Psum_" << out_lat << "_" << out_lon << "_S.dat";
    cout << "Nname: " << summed_N_name.str() << endl;
    cout << "Sname: " << summed_S_name.str() << endl;
    write_p_array(dA_N, summed_N_name.str());
    write_p_array(dA_S, summed_S_name.str());

    // Call the flux calculation:
    compFlux(dA_N, out_lat, out_lon, 0, out_dir, 
              J, flux_dist, alpha_dist);
    compFlux(dA_S, out_lat, out_lon, 1, out_dir, 
              J, flux_dist, alpha_dist);

}



