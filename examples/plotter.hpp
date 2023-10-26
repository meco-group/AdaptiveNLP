#ifndef __PLOTTER__
#define __PLOTTER__
#include <vector>
#include <optional>
#include <string>

class Plotter{
    public:
        Plotter(){};

        Plotter(bool free_time, int nx, int nu, double view_radius=0.0);

        Plotter(bool free_time, int nx, int nu, int N, int final_ind,
                std::vector<std::optional<int>>& next_ind,
                std::vector<std::optional<double>>& time_from_ind);

        void update(int N, int final_ind,
                    std::vector<std::optional<int>>& next_ind,
                    std::vector<std::optional<double>>& time_from_ind);

        // void plotTrajectory(const double* sol, std::string c="b", 
        //                     bool new_figure=true);
        // void plotTrajectory(std::vector<double>& xx,
        //                     std::vector<double>& yy, std::string c="b", 
        //                     bool new_figure=true);

        // void plotControls(const double* sol);
        // void plotControls(std::vector<std::vector<double>>& uu);
        
        void writeControlsToFile(const double* sol, int file_nb = -1);
        void writeControlsToFile(std::vector<std::vector<double>>& uu, 
                                 int file_nb = -1);

        // void plotTravelledTrajectory(std::vector<std::vector<double>> x0s,
        //                              std::string c, int frame_nb);

        // void plotAnimationFrame(int frame_nb);

        void plotAnimationFramePython(int frame_nb);

        // void plotControlAnimationFrame(int frame_nb);
        
        void addAdaptiveResults(
            std::vector<std::vector<std::vector<double>>>
                formatted_solutions_adaptive,
            std::vector<std::vector<double>> x0s_adaptive,
            std::vector<std::vector<bool>> visible_constraints_adaptive,
            std::vector<int> nb_constraints_adaptive,
            std::vector<std::vector<int>> nb_constraints_specific_adaptive,
            std::vector<double> timings_adaptive,
            std::vector<double> total_timings_adaptive);
        
        void addAdaptiveResults(
                std::vector<std::vector<std::vector<double>>>
                    formatted_solutions_adaptive,
                std::vector<std::vector<double>> x0s_adaptive,
                std::vector<std::vector<bool>> visible_constraints_adaptive,
                std::vector<int> nb_constraints_adaptive,
                std::vector<std::vector<int>> nb_constraints_specific_adaptive,
                std::vector<double> timings_adaptive,
                std::vector<double> total_timings_adaptive, 
                int scaling_x_value){
            addAdaptiveResults(formatted_solutions_adaptive, x0s_adaptive, 
                visible_constraints_adaptive, nb_constraints_adaptive, 
                nb_constraints_specific_adaptive, timings_adaptive, 
                total_timings_adaptive);
            scaling_Ns_adaptive_[scaling_Ns_adaptive_.size()-1] = 
                scaling_x_value;
        };

        void addAdaptiveIterCounts(std::vector<int>& iter_counts);

        void addRefResults(
            std::vector<std::vector<std::vector<double>>>
                formatted_solutions_ref,
            std::vector<std::vector<double>> x0s_ref,
            std::vector<std::vector<bool>> visible_constraints_ref,
            std::vector<int> nb_constraints_ref,
            std::vector<std::vector<int>> nb_constraints_specific_ref,
            std::vector<double> timings_ref,
            std::vector<double> total_timings_ref);

        void addRefResults(
                std::vector<std::vector<std::vector<double>>>
                    formatted_solutions_ref,
                std::vector<std::vector<double>> x0s_ref,
                std::vector<std::vector<bool>> visible_constraints_ref,
                std::vector<int> nb_constraints_ref,
                std::vector<std::vector<int>> nb_constraints_specific_ref,
                std::vector<double> timings_ref,
                std::vector<double> total_timings_ref, 
                int scaling_x_value){
            addRefResults(formatted_solutions_ref, x0s_ref, 
                visible_constraints_ref, nb_constraints_ref, 
                nb_constraints_specific_ref, timings_ref, total_timings_ref);
            scaling_Ns_ref_[scaling_Ns_ref_.size()-1] = scaling_x_value;
        };

        void addRefIterCounts(std::vector<int>& iter_counts);

        void addNaiveResults(
            std::vector<std::vector<std::vector<double>>>
                formatted_solutions_naive,
            std::vector<std::vector<double>> x0s_naive,
            std::vector<std::vector<bool>> visible_constraints_naive,
            std::vector<int> nb_constraints_naive,
            std::vector<std::vector<int>> nb_constraints_specific_naive,
            std::vector<double> timings_naive,
            std::vector<double> total_timings_naive);

        void addNaiveResults(
                std::vector<std::vector<std::vector<double>>>
                    formatted_solutions_naive,
                std::vector<std::vector<double>> x0s_naive,
                std::vector<std::vector<bool>> visible_constraints_naive,
                std::vector<int> nb_constraints_naive,
                std::vector<std::vector<int>> nb_constraints_specific_naive,
                std::vector<double> timings_naive,
                std::vector<double> total_timings_naive, 
                int scaling_x_value){
            addNaiveResults(formatted_solutions_naive, x0s_naive, 
                visible_constraints_naive, nb_constraints_naive, 
                nb_constraints_specific_naive, timings_naive, 
                total_timings_naive);
            scaling_Ns_naive_[scaling_Ns_naive_.size()-1] = scaling_x_value;
        };

        void addNaiveIterCounts(std::vector<int>& iter_counts);

        // void plotComputationTimes(std::vector<double>& timings);

        // void plotComputationTimes(std::vector<std::vector<double>>& timings);

        // void plotComputationTimes();

        // void plotComputationTimesFull();

        void plotComputationTimesFullFull();

        // void plotIterCounts();

        // void plotScaling();

        void writeScaling();

        void addCircle(double x, double y, double r);

        // x_min, x_max, y_min, y_max
        void addCorridor(std::vector<double>& corridor);

        void setCircleVisibility(int ind, bool visibility);

        void setCircleColor(int ind, std::string color);

        void formatSolution(std::vector<double>& sol, 
                            std::vector<std::vector<double>>& formatted_xx,
                            std::vector<std::vector<double>>& formatted_uu);
        
        void printX0s();
        
        double getMedian(std::vector<double>& v);

        void setOutputFolder(std::string output_folder){
            output_folder_ = output_folder;
        }

    private:

        void get_x(const double* sol, int k);

        void get_u(const double* sol, int k);

        // void plotCircle(int ind);

        // void plotCircle(double x, double y, double r, std::string c);

        // void plotCorridor(int ind);

        // void plotVehicle(std::vector<double>& corner_points_x, 
                        //  std::vector<double>& corner_points_y);

        std::string getMedianText(std::vector<double>& v);

        std::string output_folder_ = "../../examples/plotting_data/";

        bool free_time_;
        int nx_;
        int nu_;
        double view_radius_;
        int N_;
        int final_ind_;
        std::vector<std::optional<int>> next_ind_;
        std::vector<double> time_x_;
        std::vector<double> time_u_;

        std::vector<double> xk_;
        std::vector<double> uk_;

        std::vector<double> circle_x_;
        std::vector<double> circle_y_;
        std::vector<double> circle_r_;
        std::vector<double> visible_circles_;
        std::vector<std::string> circle_colors_;

        std::vector<std::vector<double>> corridors_;

        int plotter_count_ = 0;

        // adaptive simulation solution
        bool adaptive_solution_available_ = false;
        std::vector<std::vector<std::vector<double>>>
            formatted_solutions_adaptive_;
        std::vector<std::vector<double>> x0s_adaptive_;
        std::vector<std::vector<bool>> visible_constraints_adaptive_;
        std::vector<int> nb_constraints_adaptive_;
        std::vector<std::vector<int>> nb_constraints_specific_adaptive_;
        std::vector<double> timings_adaptive_;
        std::vector<double> total_timings_adaptive_;
        std::vector<int> iter_counts_adaptive_;
        bool adaptive_iter_counts_available_ = false;
        std::vector<double> scaling_median_total_times_adaptive_;
        std::vector<double> scaling_median_times_adaptive_;
        std::vector<double> scaling_Ns_adaptive_;

        // ref simulation solution
        bool ref_solution_available_ = false;
        std::vector<std::vector<std::vector<double>>>
            formatted_solutions_ref_;
        std::vector<std::vector<double>> x0s_ref_;
        std::vector<std::vector<bool>> visible_constraints_ref_;
        std::vector<int> nb_constraints_ref_;
        std::vector<std::vector<int>> nb_constraints_specific_ref_;
        std::vector<double> timings_ref_;
        std::vector<double> total_timings_ref_;
        std::vector<int> iter_counts_ref_;
        bool ref_iter_counts_available_ = false;
        std::vector<double> scaling_median_total_times_ref_;
        std::vector<double> scaling_median_times_ref_;
        std::vector<double> scaling_Ns_ref_;

        // naive simulation solution
        bool naive_solution_available_ = false;
        std::vector<std::vector<std::vector<double>>>
            formatted_solutions_naive_;
        std::vector<std::vector<double>> x0s_naive_;
        std::vector<std::vector<bool>> visible_constraints_naive_;
        std::vector<int> nb_constraints_naive_;
        std::vector<std::vector<int>> nb_constraints_specific_naive_;
        std::vector<double> timings_naive_;
        std::vector<double> total_timings_naive_;
        std::vector<int> iter_counts_naive_;
        bool naive_iter_counts_available_ = false;
        std::vector<double> scaling_median_total_times_naive_;
        std::vector<double> scaling_median_times_naive_;
        std::vector<double> scaling_Ns_naive_;

        int dpi_ = 100;
        std::string file_format_ = ".png";
};

#endif