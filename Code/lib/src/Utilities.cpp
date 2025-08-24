#include "Utilities.h"
#include <chrono>
#include <cassert>

Rnd rnd;



////////////////// always implemented: //////////////////

// finds x for which func(x) is lowerr than tol, !!! make sure func crosses zero in range x_left xright at least and only one time !!!
double bisection(double x_left, double x_right, std::function<double(double)>& func){
    constexpr int N_max = 64;
    constexpr double tol = 1e-7;
    
    double y_left = func(x_left);
    double y_right = func(x_right);
    
    if(y_left < 0 && y_right > 0)
    {
        for(int i = 0; i < N_max; i++){
            double x_new = (x_left + x_right)*0.5;
            double y_new = func(x_new);

            dmsg( "y_new: " << y_new << "\n");
            if( std::abs(y_new) < tol )
            {
                dmsg("returning normally case 1\n");
                return x_new;
            }
            
            if(y_new < 0)
            {
                x_left = x_new;
            }
            else
            {
                x_right = x_new;
            }
        }
    }
    else if(y_right < 0 && y_left > 0)
    {
        for(int i = 0; i < N_max; i++){
            double x_new = (x_left + x_right)*0.5;
            double y_new = func(x_new);

            dmsg( "y_new: " << y_new << "\n");
            if( std::abs(y_new) < tol )
            {
                dmsg("returning normally case 2\n");
                return x_new;
            }
            
            if(y_new > 0)
            {
                x_left = x_new;
            }
            else
            {
                x_right = x_new;
            }
        }
    }
    else{
        throw std::invalid_argument("Error in bisection! func(x_left) == func(x_right) = 0 or both have same sign: " +str(y_left) + " | " + str(y_right) + "\n");
    }

    dmsg("returning case 3\n");
    return (x_left + x_right) * 0.5;

};

double newton(int n, int k){
    if(n < 0 || k < 0) return 0;
    if(n < 160 && (k < 13 || n-k < 13)) return newton_size_t(n,k);
    if(n < 271 && (k < 11 || n-k < 11)) return newton_size_t(n,k);
    if(n==63 && k>28 && k<35) return newton_double(n,k);
    if(n==64 && k>26 && k<38) return newton_double(n,k);
    if(n==65 && k>25 && k<40) return newton_double(n,k);
    if(n==66 && k>24 && k<42) return newton_double(n,k);
    if(n==67 && k>23 && k<44) return newton_double(n,k);
    if(n==68 && k>23 && k<45) return newton_double(n,k);
    if(n==69 && k>22 && k<47) return newton_double(n,k);
    if(n > 1020){
        double newton_value = newton_double(n,k);
        assert(!std::isinf(newton_value));
        return newton_value;
    }   
    if(n>69) return newton_double(n,k);
    return newton_size_t(n,k); 
}

double newton_double(int n, int k){ // wont overflow until n > 1020 and k > 495
    if(k>n) return 0;
    if(k>n-k) k = n-k;
    double out = 1;

    for(int i = 1; i <=k; i++){
        out *= n;
        out /= i;
        n--;
    }
    return out;
};

size_t newton_size_t(int n, int k){ // will overflow fast but if not than 100% accurate
    if(k>n) return 0;
    if(k>n-k) k = n-k;
    size_t out = 1;

    for(int i = 1; i <= k; i++){
        out *= n;
        out /= i;
        n--;
    }
    return out;
};

double time(){ // returns time in seconds
    static auto time_start = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::high_resolution_clock::now();
    return (time - time_start).count()*1e-9;
};

Rnd::Rnd(): mt_gen{std::random_device()()}, rnd_dist{0.0, 1.0} { //[0, 1.0)
};
Rnd::Rnd(unsigned seed): mt_gen{seed}, rnd_dist{0.0, 1.0} { //[0, 1.0)
};

double Rnd::operator()(){
    return rnd_dist(mt_gen);
};
double Rnd::operator()(const double& min, const double& max){
    return min + rnd_dist(mt_gen)*(max - min); 
};

bool save(double data, const std::string& path, const char mode){
    return save_scalar(data, path, mode);
}
bool save(size_t data, const std::string& path, const char mode){
    return save_scalar(data, path, mode);
}
bool save(int data, const std::string& path, const char mode){
    return save_scalar(data, path, mode);
}
void reset_color(){
    std::cout << "\033[0m ";
}
void reset_color(std::stringstream& ss){
    ss << "\033[0m ";
}


size_t permute_transpose_inplace(size_t a, size_t new_first_dim, size_t new_second_dim){
    size_t MN_1 = new_first_dim*new_second_dim-1; 
    if(a == MN_1) return MN_1;
    return (new_first_dim*a)%MN_1;  
};

bool is_minimal_in_cycle(size_t a, size_t new_first_dim, size_t new_second_dim){
    // size_t MN_1 = new_first_dim*new_second_dim-1; 
    size_t current_value = permute_transpose_inplace(a, new_first_dim, new_second_dim);

    while(current_value!=a){
        if(current_value<a) return false;
        current_value = permute_transpose_inplace(current_value, new_first_dim, new_second_dim);
    }
    return true;
};

size_t inverse_permute_transpose_inplace(size_t a, size_t new_first_dim, size_t new_second_dim){
    size_t MN_1 = new_first_dim*new_second_dim-1; 
    if(a == MN_1) return MN_1;
    return (new_second_dim*a)%MN_1;  
};

// a b c
// d e f
// g h i
// j k l

// a b c d e f g h  i  j k l
// 0 4 8 1 5 9 2 6  10 3 7 11 inverse permute orig -> dest na ten index idzie ten element
// 0 3 6 9 1 4 7 10 2  5 8 11 permute         dest <- orig na to miejsce przychodzi taki index
// a d g j b e h k  c  f i l

// 0 1 2 3 4 5 6 7  8  9 10 11
// a d g
// j b e
// h k c
// f i l

// a d g j
// b e h k
// c f i l

// a d g j
// b e h k
// c f i l
