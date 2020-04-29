////////////////////////////////////////////////////////////////////////////////
// This program generates 10 stock price paths using a GARCH(1,1)
//    stochastic volatility model. Results are printed to the screen. The
//    current (stochastic) and long-run volatility are 30% annually.
////////////////////////////////////////////////////////////////////////////////

#include "Functions.h"


// Used for calculating the max payoff of a simple vanilla call
double Max(double& val) {
    if (val > 0.0) { return val; }
    return 0.0;
}

// This is a function for plotting x and y data in TeX
void Plot(double x[], double y[]) {
    FILE* fp;
    fp = fopen("Plot.tex", "w");

    // construct the plot
    fprintf(fp, "\\documentclass{standalone}\n");
    fprintf(fp, "\\usepackage{pgfplots}\n");
    fprintf(fp, "\\begin{document}\n");
    fprintf(fp, "\\begin{tikzpicture}\n");
    // add a plot 
    fprintf(fp, "\\begin{axis}[\n");
    fprintf(fp, "title={Volatility Skew},\n");
    fprintf(fp, "xlabel={Strikes},\n");
    fprintf(fp, "ylabel={Implied Volatility (IV)},\n");
    fprintf(fp, "xmin=%8.4f, xmax=%8.4f,\n", x[0] * .80, x[10] * 1.20);
    fprintf(fp, "ymin=%8.4f, ymax=%8.4f,\n", y[0] * .80, y[10] * 1.20);
    fprintf(fp, "ymajorgrids=true,\n");
    fprintf(fp, "grid style=dashed,\n");
    fprintf(fp, "]\n");
    fprintf(fp, "\\addplot ");

    fprintf(fp, "coordinates {\n");
    for (int i = 0; i < 11; ++i) {
        fprintf(fp, "(%8.4f, %8.4f)\n", x[i], y[i]);
    }
    fprintf(fp, "};\n");
    // complete the construction 
    fprintf(fp, "\\end{axis}\n");
    fprintf(fp, "\\end{tikzpicture}\n");
    fprintf(fp, "\\end{document}\n");
    fclose(fp);

    return;
}

int main() {

    int i, j, n, N, sims;
    double s_start, s, T, r, sigma, sigma2, sigma2_start, mu, alpha, beta, r_annual,
        gamma, dt, N01, R, U, V2, K, val1, val2, avg1;

    // Define the strikes, as well as declare the averages (Problem 2), and volatility (Problem 3) arrays. 
    double strikes[11] = { 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0 };
    double averages[11] = { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 };
    double vol[11] = { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 };

    // alpha, beta, and gamma parameters as computed in the XOM example.
    alpha = 0.045; beta = 0.083; gamma = 0.872;

    // The average terminal stock price we would like to calculate for the problem 1 
    avg1 = 0.0;

    // Number of simulations to run
    sims = 10000000;

    // Time to expiration
    T = 0.5;

    // Risk-free interest rate.
    r_annual = 0.05; // annualized rate
    r = r_annual;    // rate that will be converted into daily rate

    // Number of days to expiration (252 days per year), and one day in years.
    N = 252 * T;
    dt = 1.0 / 252.0;

    // Convert "r" to a daily risk-free rate.
    r *= dt;

    // Current stock price.
    s_start = 100.0;

    // Current stock price variance.
    sigma2_start = 0.35 * 0.35;

    // Current daily stock price variance.
    sigma2_start *= dt;

    // Annual long-term variance.
    V2 = 0.30 * 0.30;

    // Daily long-term variance.
    V2 *= dt;

    // Seed the RNG.
    MTUniform(1);

    // Generate Monte Carlo paths for Stochastic volatility
    for (n = 1; n <= sims; ++n) {

        // Initialize the stock price and the volatility to their starting (time 0) values.
        s = s_start;
        sigma2 = sigma2_start;

        // Now march forward in time day by day, generating stock price path
        for (i = 1; i <= N; ++i) {

            // Compute the drift parameter that corresponds to the volatility for this period.
            mu = r - (0.5 * sigma2);

            // Compute the stock price at day "i"...

            // First get a standard normal N01.
            U = MTUniform(0);
            N01 = PsiInv(U);

            // Apply current volatility.
            sigma = sqrt(sigma2);
            R = sigma * N01;

            // Update the stock price.
            s *= exp(mu + R);

            // Update the stochastic volatility according to the GARCH(1,1) model.
            sigma2 = alpha * V2 + beta * R * R + gamma * sigma2;
        }

        // Problem 1 //
        // Since K = 0, we simply need to find the average terminal stock price.
        // Note: "avg1" isn't the actual average yet, we still need to divide it by the number of simulations
        avg1 += exp(-r_annual * T) * s;

        // Problem 2 //
        for (j = 0; j < 11; ++j) {
            // Iterate through each of the strike values
            // Subtract the strike from terminal stock price, then discount to T=0.
            // The Max function takes the maximum between 0, and the number provided
            // Note: "averages" isn't the actual average yet, we still need to divide it by the number of simulations
            val2 = (s - strikes[j]) * exp(-r_annual * T);
            averages[j] += Max(val2);
        }

        if (n % 100000 == 0) {
            printf("%3d percent of the simulations have been completed...\n", 100 * n / sims);
        }

    }

    // Outputs value for problem 1
    printf("\nProblem 1:\n");
    printf("Our call option with strike K = 0 is valued at %8.4f\n", avg1 / sims);
    printf("--------------------------------------------------------------\n");

    // Outputs value for problem 2

    printf("Problem 2:\n");
    for (int t = 0; t < 11; ++t) {
        printf("A call option with strike K = %8.2f is valued at %8.4f\n", strikes[t], averages[t] / sims);
    }
    printf("--------------------------------------------------------------\n");

    // Outputs value for problem 3

    // define x, y arrays
    printf("Problem 3:\n");
    for (int t = 0; t < 11; ++t) {
        K = strikes[t];
        // Function inputs: time to expiry, starting stock price, strike price, annual risk-free rate, option price
        vol[t] = 100.0 * ImpliedVol(T, s_start, K, r_annual, averages[t] / sims);

        printf("Our call option with strike K = %8.2f has implied volatility at %8.4f%\n", K, vol);
    }
    printf("--------------------------------------------------------------\n");

    // Outputs a TeX file for plotting
    Plot(strikes, vol);

    Exit();
}


