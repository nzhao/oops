#include <vector>

class Spin
{
public:
    Spin(std::vector<double> coordinate)
    {
    }
    std::vector<double> get_coordinate();
    void set_coordinate(std::vector<double>);
    
private:
    std::vector<double> coordinate;
};
