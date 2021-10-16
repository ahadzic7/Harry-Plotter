//Faculty of Electrical Engineering Sarajevo
//My simple plotter
//Armin Hadzic

#include "mbed.h"
#include "stm32f413h_discovery_ts.h"
#include "stm32f413h_discovery_lcd.h"
#include "cmath"
#include <vector>
#include <functional>
#include <algorithm>
typedef std::pair<double, double> Par;
typedef std::vector<Par> VectorPar;
typedef std::vector<double> Vector;

bool jednakost(double a, double b) { return std::fabs(a - b) < 1e-6; }

bool fun_sort(const Par &a, const Par &b) {
    // if (jednakost(a.first, b.first))
    //     throw std::domain_error("Invalid data set");
    return a.first < b.first;
}

bool fun_criteria(const Par &a, const double tacka) { return a.first < tacka; }


class AbstractInterpolator {
    protected:
        typedef std::pair<Par, Par> ParParova;

        int Locate(double x) const {
            if (x < tacke[0].first || jednakost(x, tacke[0].first))
                return 0;
            if (x > tacke[tacke.size() - 1].first || jednakost(tacke[tacke.size() - 1].first, 0))
                return tacke.size();

            if(provjeraKesha(x, keshIndex))
                return keshIndex + 1;
            if(keshIndex > 0 && provjeraKesha(x, keshIndex - 1))
                return keshIndex;
            if(keshIndex < tacke.size() - 1  && provjeraKesha(x, keshIndex + 1))
                return keshIndex + 2;

            keshIndex = std::lower_bound(tacke.begin(), tacke.end(), x,fun_criteria) - tacke.begin() - 1;

            return keshIndex + 1;
        }

        VectorPar getTacke() { return tacke; }
        Par getTacka(unsigned int index) const { return tacke[index]; }
        unsigned int getTackeSize() const { return tacke.size(); }
        // static bool jednakost(double a, double b) {
        //     const double EPS (0.0000001);
        //     return std::fabs(a - b) < EPS;
        // }
        void addTacka(const Par novaTacka) {  tacke.push_back(novaTacka); }

    private:
        VectorPar tacke;
        mutable int keshIndex;

        bool provjeraKesha (double x, unsigned int index) const { return x > tacke[index].first && x < tacke[index + 1].first; }

    public:
        virtual double operator()(double x) const = 0;

        explicit
        AbstractInterpolator(const VectorPar &data) : tacke(data), keshIndex(data.size() / 2) { std::sort(tacke.begin(), tacke.end(), fun_sort); }
};

class LinearInterpolator : public AbstractInterpolator {
    private:
        ParParova interval(double x) const {
            int indexTacke(Locate(x));
            if (indexTacke == 0) 
                return std::make_pair(getTacka(0), getTacka(1));
            
            if (indexTacke == getTackeSize()) 
                return std::make_pair(getTacka(indexTacke - 2), getTacka(indexTacke - 1));
            
            return std::make_pair(getTacka(indexTacke - 1), getTacka(indexTacke));
        }

    public:
        LinearInterpolator(const VectorPar & data) : AbstractInterpolator(data) {}
        double operator()(double x) const override {
            ParParova par(interval(x));
            return par.first.second + (par.second.second - par.first.second) * (x - par.first.first) / (par.second.first - par.first.first);
        }
};

class TrigonometricInterpolator : public AbstractInterpolator {
    private:
        int M;
    public:
        TrigonometricInterpolator(const VectorPar &data) : AbstractInterpolator(data) {
            // if (!jednakost(getTacka(0).second, getTacka(getTackeSize() - 1).second))
            //     throw std::domain_error("Function is not periodic");
            if(getTackeSize() % 2 == 0)
                M = (getTackeSize() - 2) / 2;
            else
                M = (getTackeSize() - 1) / 2;
        }

        double operator()(double x) const override {
            const double W (2 * M_PI / (getTacka(getTackeSize() - 1)).first - getTacka(0).first);
            double suma(0);
            if (getTackeSize() % 2 == 0) {
                for (int k = 0; k < 2 * M + 1; k++){
                    double proizvod(1);
                    for(int j = 0; j < 2 * M + 1; j++) {
                        if (j != k) {
                            proizvod *= std::sin(W / 2 * (x - getTacka(j).first)) / std::sin(W / 2 * (getTacka(k).first - getTacka(j).first));
                        }
                    }
                    suma += getTacka(k).second * proizvod;
                }
            }
            else {
                for(int k = 0; k < 2 * M; k++){
                    double alfa(0);
                    for(int j = 0; j < 2 * M; j++) {
                        if(j != k)
                            alfa += getTacka(j).first;
                    }
                    alfa = -alfa;

                    double mnozioc(getTacka(k).second * std::sin(W / 2 * (x - alfa)) / std::sin(W / 2 * (getTacka(k).first - alfa)));

                    double proizvod(1);
                    for(int j = 0; j < 2 * M; j++) {
                        if (j != k) 
                            proizvod *= std::sin(W / 2 * (x - getTacka(j).first)) / std::sin(W / 2 * (getTacka(k).first - getTacka(j).first));
                    }
                    suma += mnozioc * proizvod;
                }
            }
            return suma;
        }
};

class SplineInterpolator : public AbstractInterpolator {
    private:
        Vector r, s, q;
    public:
        explicit
        SplineInterpolator(const VectorPar &data) : AbstractInterpolator(data), r(data.size()), s(data.size() - 1), q(data.size() - 1) {
            const int VELICINA(getTackeSize());

            for(int i = 1; i < VELICINA - 1; i++) {
                q[i] = 2 * (getTacka(i + 1).first - getTacka(i - 1).first);
                r[i] = 3 * ( (getTacka(i + 1).second - getTacka(i).second) / (getTacka(i + 1).first - getTacka(i).first)
                           - (getTacka(i).second - getTacka(i - 1).second) / (getTacka(i).first - getTacka(i - 1).first) );
            }

            for(int i = 1; i < VELICINA - 2; i++) {
                double mi ((getTacka(i + 1).first - getTacka(i).first) / q[i]);
                q[i + 1] -= mi * (getTacka(i + 1).first - getTacka(i).first);
                r[i + 1] -= mi * r[i];
            }

            r[VELICINA - 2] /= q[VELICINA - 2];

            for(int i = VELICINA - 3; i > 0; i--) 
                r[i] = (r[i] - (getTacka(i + 1).first - getTacka(i).first) * r[i + 1]) / q[i];

            for(int i = 0; i < VELICINA - 1; i++) {
                double dx(getTacka(i + 1).first - getTacka(i).first);
                s[i] = (r[i + 1] - r[i]) / (3 * dx);
                q[i] = (getTacka(i + 1).second - getTacka(i).second) / dx - dx * (r[i + 1] + 2 * r[i]) / 3;
            }
        }

        double operator()(double x) const override {
            int index(Locate(x));
            if (index != 0)
                index--;
            if(index == getTackeSize() - 1)
                index--;
            double t (x - getTacka(index).first);
            return getTacka(index).second + t * (q[index] + t * (r[index] + t * s[index]));
        }
};

class BarycentricInterpolator : public AbstractInterpolator {
    private:
        Vector w;
    public:
        BarycentricInterpolator(const VectorPar &data, int order) : AbstractInterpolator(data), w(data.size()) {
            const int N(getTackeSize());
            if (order < 0 || order > N)
                throw std::domain_error("Invalid order");

            for(int i = 0; i < N; i++) {
                w[i] = 0;
                double p;
                for(int k = std::max(1, i - order); k <= std::min(i, N - order) + 1; k++) {
                    p = 1;
                    for(int j = k; j < k + order; j++){
                        if (j != i) {
                            p = p / (getTacka(i).first - getTacka(j).first);
                        }
                    }
                    if (k % 2 == 0) {
                        p = -p;
                    }
                }
                w[i] += p;
            }
        }

        double operator()(double x) const override {
            const int N(getTackeSize());
            double p = 0;
            double q = 0;

            for(int i = 0; i < N; i++) {
                if (jednakost(x, getTacka(i).first))
                    return getTacka(i).second;

                double u (w[i] / (x - getTacka(i).first));
                p += u * getTacka(i).second;
                q += u;
            }
            return p / q;
        }

        Vector GetWeights() { return w; }
};

class PiecewisePolynomialInterpolator : public AbstractInterpolator {
    private:
        int red;

        double lagrangePolinom(const VectorPar &vek, double x) const {
            double s(0);
            const int VELICINA(vek.size());

            for(int i = 0; i < VELICINA; i++) {
                double p(vek[i].second);
                for(int j = 0; j < VELICINA; j++) {
                    if (i != j)
                        p *= (x - vek[j].first) / (vek[i].first - vek[j].first);
                }
                s += p;
            }
            return s;

        }
        double lagrangePolinom2(int left, int right, double x) const {
            double s(0);
            const int VELICINA(right - (left - 1));
            const int pom2(left);

            for(int i = 0; i < VELICINA; i++) {
                double p(getTacka(left).second);
                int pom(pom2);
                for(int j = 0; j < VELICINA; j++) {
                    if (i != j)
                        p *= (x - getTacka(pom).first) / (getTacka(left).first - getTacka(pom).first);
                    pom++;
                }
                left++;
                s += p;
            }
            return s;

        }

    public:
        PiecewisePolynomialInterpolator(const VectorPar &data, const int order) : AbstractInterpolator(data) {
            if(order <= 1 || order > getTackeSize())
                throw std::domain_error("Invalid order");
            red = order;
        }


        double operator()(double x) const override {
            int index(Locate(x) -1);
            const int VELICINA(getTackeSize());
            int lijevaGranica, desnaGranica;

            if (index == VELICINA - 1)
                index--;

            if(red % 2 == 0) {
                lijevaGranica = index - red / 2;
                desnaGranica = index + red / 2;
            } else {
                lijevaGranica = index - (red - 1) / 2;
                desnaGranica = index + (red + 1) / 2;
            }
            if(lijevaGranica < 1) {
                lijevaGranica = 0;
                desnaGranica = red;
            }
            else if (desnaGranica >= VELICINA) {
                lijevaGranica = desnaGranica - (red + 1);
                desnaGranica = VELICINA - 1;
            }
            const int VEL (desnaGranica - (lijevaGranica - 1));
            VectorPar y(VEL);

            for(int i = 0; i < VEL; i++) 
                y[i] = getTacka(lijevaGranica);

            return lagrangePolinom2(lijevaGranica, desnaGranica, x);
        }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Dot {
    double x;
    double y;
    Dot(double x, double y) { this->x = x; this->y = y;}
};

TS_StateTypeDef TS_State = { 0 };
AnalogIn A(p15);
AnalogIn B(p16);
AnalogIn C(p17);
AnalogIn S(p18);
AnalogIn X(p19);
double a, b, c, scale;

bool scale_change() {
    if(jednakost(scale, S.read())) return false;
    scale = S.read();
    return true;
}
bool parametars_change() {
    if(jednakost(a, A.read()) && jednakost(b, B.read()) && jednakost(c, C.read())) return false;
    a = A.read();
    b = B.read();
    c = C.read();
    return true;
}

AbstractInterpolator *interpolator;

LinearInterpolator lin(VectorPar(3));
PiecewisePolynomialInterpolator pin(VectorPar(3), 3);
TrigonometricInterpolator tin(VectorPar(3));
BarycentricInterpolator bin(VectorPar(3), 3);
SplineInterpolator spin(VectorPar(3));

double get_interpolate(double x) { return interpolator->operator()(x); }
double get_sin(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = M_PI * (c - 0.5);
    return a_ * sin(b_ * x + c_); 
}
double get_cos(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = M_PI * (c - 0.5);
    return a_ * cos(b_ * x + c_); 
}
double get_sinc(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = M_PI * (c - 0.5);

    if(x<-.001 || x>.001)
        return a_ * sin(b_ * x + c_) / x;
    return a_ * sin(b_ * 0.002 + c_) / 0.002;
}
double get_saw(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = M_PI * (c - 0.5);

    return a_ - (2 * a_) / M_PI * atan(1. / tan((M_PI * x) / b_)); 
}
double get_triangle(double x) {
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = M_PI * (c - 0.5);
    return a * acos(cos(b * x + c));
}
double get_squares(double x) { 
    double y(get_sin(x));
    if(y < 0) return -1;
    return 1;
}
double get_exponential(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b);
    double c_ = 5 * (c - 0.5);

    return a_ * std::pow(b_, x + c_); 
}
double get_logarithmic(double x) {
    if(x < 0) return -1./0;
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b + 0.5);
    double c_ = 5 * (c - 0.5);

    return a_ * log(b_ * x) + c_; 
}
double get_linearf(double x) { 
    double a_ = (a - 0.5) * 3.14;
    double b_ = (b - 0.5) * 5;
    return a_ * x + b_; 
}
double get_squaref(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = 5 * (c - 0.5);
    return a_ * x * x + b_ * x + c_; 
}
double get_cubef(double x) { 
    double a_ = 5 * (a - 0.5);
    double b_ = 5 * (b - 0.5);
    double c_ = 5 * (c - 0.5);
    return a * x * x * x + b * x * x + c * x;
}
double get_inversef(double x) { return 1 / x; }

struct Pixel {
    unsigned int x, y;
    Pixel(int x, int y) { this->x = x; this->y = y; }
    bool operator == (const Pixel &p2) { return x == p2.x && y == p2.y; }
    bool operator < (const Pixel &p2) {
        if(x == p2.x) return y < p2.y;
        return x < p2.x;
    }
};
//////////////////////////////////////////////////////////////////////////////
class Plotter {
    private:
        static const unsigned int DIM = 240;
        static const unsigned int MAX_SCALE = 18;
        static const unsigned int DOTR = 3;
        int widthX;
        Dot upRight, downLeft;
        std::function<double(double)> function;
        std::vector<Pixel> pixels;
        VectorPar interpolateDots;

        unsigned int bgColour, fColour, cooColour, xColour, gColour, tColour, dColour;
        unsigned int interCounter;

        unsigned short state, pause;

        void line(Pixel p1, Pixel p2) { BSP_LCD_DrawLine (p1.x, p1.y, p2.x, p2.y); }
        void hline(Pixel t, int length) { BSP_LCD_DrawHLine (t.x, t.y, length); }
        void vline(Pixel t, int length) { BSP_LCD_DrawVLine (t.x, t.y, length); }

        void coordinate_system() {
            BSP_LCD_SetTextColor(cooColour);
            vline(Pixel(DIM / 2, 0), DIM);
            hline(Pixel(0, DIM / 2), DIM);
            BSP_LCD_SetTextColor(fColour);
        }

        void grid() {
            for(int i = 0; i < DIM; i += DIM / widthX) {
                for(int j = 0; j < DIM; j += DIM / widthX) 
                    BSP_LCD_DrawPixel(i, j, gColour);
            }
        }

        int rtpX(double x) { //converts real to pixels on X axis
            if(x < downLeft.x || x > upRight.x) 
                return -1;
            return DIM * (x + upRight.x) / widthX;
        }
        
        int rtpY(double y) { //converts real to pixels on Y axis
            if(y < downLeft.y || y > upRight.y) 
                return -1;
            return DIM * (1 - (y + upRight.x) / widthX);
        }
        //inverse
        double ptrX(unsigned int pix) { return static_cast<double>((pix * widthX)) / DIM - upRight.x; }

        double ptrY(unsigned int pix) { return static_cast<double>((DIM - pix) * widthX) / DIM - upRight.x; }

        int toScale() { return S.read() * 18 + 2; }

        //calculate function values
        std::vector<Pixel> values() {
            double x (downLeft.x);
            std::vector<Pixel> pixels;
            for(int i = 0; i < DIM; i++) {
                double y = function(x);
                pixels.push_back(Pixel(i, rtpY(y)));
                x = ptrX(i);
            }
            return pixels;
        }

        //displays key parametars to the user
        void display_value(double x, double y) {
            BSP_LCD_SetTextColor(LCD_COLOR_WHITE);

            BSP_LCD_FillRect(0, 224, 240, 240);
            BSP_LCD_FillRect(0, 0, 240, 16);

            BSP_LCD_SetTextColor(fColour);
            BSP_LCD_SetFont(&Font16);
            char value[20] = {'\0'};
            if(x < 10 && y < 10)
                sprintf(value, "x = %.2f y = %.2f", x, y);
            if(x < 100 && y < 100)
                sprintf(value, "x = %.1f y = %.1f", x, y);
            else 
                sprintf(value, "x = %.0f y = %.0f", x, y);

            char params[20] = {'\0'};
            sprintf(params, "%dx%d(%.2f,%.2f,%.2f)", int(upRight.x), int(upRight.x), a, b, c);

            hline(Pixel(0, 16), DIM);
            hline(Pixel(0, 224), DIM);

            BSP_LCD_SetTextColor(tColour);
            BSP_LCD_SetBackColor(LCD_COLOR_WHITE);
            BSP_LCD_DisplayStringAt(0, 0, (uint8_t *)params, CENTER_MODE);
            BSP_LCD_DisplayStringAt(0, 225, (uint8_t *)value, CENTER_MODE);
            BSP_LCD_SetTextColor(fColour);
        }
        //display for interpolate dots input
        void interpolate_display() {

            BSP_LCD_Clear(bgColour);

            grid();
            coordinate_system();
            
            BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
            BSP_LCD_SetBackColor(LCD_COLOR_WHITE);


            BSP_LCD_FillRect(0, 224, 240, 240);
            BSP_LCD_FillRect(0, 0, 240, 16);

            BSP_LCD_SetTextColor(fColour);
            BSP_LCD_SetFont(&Font16);
            char message[] = "Unesite tacke";
            char message2[] = "Minimalno dvije tacke";
          

            hline(Pixel(0, 16), DIM);
            hline(Pixel(0, 224), DIM);

            BSP_LCD_SetTextColor(tColour);
            BSP_LCD_DisplayStringAt(0, 0, (uint8_t *)message, CENTER_MODE);
            BSP_LCD_DisplayStringAt(0, 225, (uint8_t *)message2, CENTER_MODE);

            for(Par par: interpolateDots) {
                Pixel pix(rtpX(par.first), rtpY(par.second));
                BSP_LCD_SetTextColor(LCD_COLOR_RED);
                BSP_LCD_FillCircle(pix.x, pix.y, 3);
            }

            BSP_LCD_SetTextColor(fColour);
        }
        //draws the actual function
        void draw(bool interpolate) {
            if(!(state == 3 || state == 6)) return;
            bool live(false);
            //get new values if somthing is changed
            if(scale_change()) {
                set_scale(toScale());
                pixels = values();
                live = true;
            }
            else if(parametars_change()) {
                pixels = values();
                live = true;
            }
            else if(interpolate && interCounter == 0) {
                pixels = values();
                live = true;
                interCounter++;
            }
            BSP_LCD_SetTextColor(dColour);
            for(Par p: interpolateDots)
                BSP_LCD_FillCircle(rtpX(p.first), rtpY(p.second), DOTR); 
        
            BSP_LCD_SetTextColor(fColour);
            for(int i = 0; i < pixels.size() - 2; i++) {
                if(!(state == 3 || state == 6)) return;
                line(pixels[i], pixels[i + 1]);
                if(live && get_pause() > 0) wait_ms(get_pause());
            }
        }


    public:

        Plotter(Dot uR, Dot dL) : widthX(uR.x - dL.x), upRight(uR), downLeft(dL), state(0), pause(0), interCounter(0) {
            bgColour = LCD_COLOR_WHITE;
            fColour = LCD_COLOR_BLACK;
            cooColour = LCD_COLOR_RED;
            xColour = LCD_COLOR_GREEN;
            gColour = LCD_COLOR_GRAY;
            tColour = LCD_COLOR_BLACK;
            dColour = LCD_COLOR_RED;
        }

        void screen_setup() {
            if(!(state == 3 || state == 6)) return;
            BSP_LCD_Init();
            if (BSP_TS_Init(BSP_LCD_GetXSize(), BSP_LCD_GetYSize()) == TS_ERROR) 
                printf("BSP_TS_Init error\n");
            
            
            BSP_LCD_Clear(bgColour);
            BSP_LCD_SetTextColor(fColour);
            grid();
            coordinate_system();
        }
        
        void set_scale(int scale) { upRight = Dot(scale, scale); downLeft = Dot(-scale, -scale);  widthX = 2 * scale; }

        void set_function(const std::function<double(double)> &function) { this->function = function; interCounter = 0; }
        
        void set_pixels(const std::vector<Pixel> &pixels) { this->pixels = pixels; }

        void set_state(unsigned short state) { this->state = state; }

        unsigned short get_state() { return state; }

        void set_background(unsigned int colour) { 
            switch(colour) {
                case LCD_COLOR_BLACK: {
                    fColour = LCD_COLOR_WHITE;
                    gColour = LCD_COLOR_GRAY;
                    break;
                }
                case LCD_COLOR_RED: {
                    cooColour = LCD_COLOR_BLUE;
                    gColour = LCD_COLOR_WHITE;
                    dColour = LCD_COLOR_BLUE;
                    break;
                }
                case LCD_COLOR_GREEN: {
                    xColour = LCD_COLOR_BLUE;
                    gColour = LCD_COLOR_RED;
                    break;
                }
                case LCD_COLOR_GRAY: {
                    gColour = LCD_COLOR_BLACK;
                    break;
                }
                case LCD_COLOR_WHITE: {
                    bgColour = LCD_COLOR_WHITE;
                    fColour = LCD_COLOR_BLACK;
                    cooColour = LCD_COLOR_RED;
                    xColour = LCD_COLOR_GREEN;
                    gColour = LCD_COLOR_GRAY;
                    tColour = LCD_COLOR_BLACK;
                    dColour = LCD_COLOR_RED;
                    return;
                } 
            }
            bgColour = colour; 
        }

        VectorPar get_interpolate_dots() { return interpolateDots; }

        void clear_interpolate_dots() { interpolateDots.clear(); }

        void set_pause(unsigned short pause) { this->pause = pause; }
        
        unsigned short get_pause() { return pause; }
        
        double get_value(double potRead, bool interpolate = false) {
            if(!(state == 3 || state == 6)) return 0;
            if(potRead < downLeft.x || potRead > upRight.x) 
                return 0;    
            
            double x (potRead * widthX - upRight.x);
            double y (function(x));
            
            BSP_LCD_Clear(bgColour);
            coordinate_system();
            draw(interpolate);
            
            BSP_LCD_SetTextColor(xColour);
            vline(Pixel(rtpX(x), 0), DIM);
            grid();
            display_value(x, y);
            return y;
        }
        //reads touch from the screen and stores interpolate dots
        void inputInterpolateDots(bool del = true) {
            interpolate_display();
            if(del) 
                interpolateDots.clear();            
            unsigned short x(-1), y(-1);
            while(state == 4) {
                BSP_TS_GetState(&TS_State);
                if(TS_State.touchDetected) {
                    unsigned short x1 = TS_State.touchX[0];
                    unsigned short y1 = TS_State.touchY[0];
                    if(x == x1 && y == y1) continue;
                    x = x1;
                    y = y1;
                    Pixel pix(x1, y1);
                    BSP_LCD_SetTextColor(LCD_COLOR_RED);
                    BSP_LCD_FillCircle(pix.x, pix.y, 3);
                    interpolateDots.push_back(std::make_pair(ptrX(pix.x), ptrY(pix.y)));
                }
                
                wait_ms(100);
            }
        }
};

int w = 5;
Dot uR(w, w), dL(-w, -w);
Plotter p(uR, dL);
/////////////////////////////////////////// INTERFACE btn_functions ////////////////////////////////////
InterruptIn power_btn(p5);
InterruptIn confirm_btn(p6);
InterruptIn back_btn(p7);

/////////////////////////////////////VIEW FUNCTIONS

// main meni

void display_title() {
    BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
    BSP_LCD_SetFont(&Font24);
    BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
    char niz[15] = "HARRY PLOTTER";
    BSP_LCD_DisplayStringAt(0, 25, (uint8_t*)niz, CENTER_MODE);

    BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
    BSP_LCD_SetFont(&Font16);
    BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
    char niz1[15] = "Dobro dosli!";
    BSP_LCD_DisplayStringAt(5, 85, (uint8_t*)niz1, CENTER_MODE);
}
void btn_settings(){
    BSP_LCD_SetTextColor(LCD_COLOR_GRAY);
    BSP_LCD_FillRect (200, 75, 30, 30);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    BSP_LCD_SetFont(&Font24);
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    char niz[2] = "#";
    BSP_LCD_DisplayStringAt(98, 80, (uint8_t*)niz, CENTER_MODE);
}
void btn_interpolation(){
    BSP_LCD_SetTextColor(LCD_COLOR_GRAY);
    BSP_LCD_FillRect ( 9, 125, 110, 105);
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    char niz[14] = "Interpolacija";
    BSP_LCD_SetFont(&Font12);
    BSP_LCD_DisplayStringAt(18, 160, (uint8_t*)niz, LEFT_MODE);
    char niz1[9] = "funkcije";
    BSP_LCD_DisplayStringAt(35, 190, (uint8_t*)niz1, LEFT_MODE);
}
void btn_functions(){
    BSP_LCD_SetTextColor(LCD_COLOR_GRAY);
    BSP_LCD_FillRect (120, 125, 110, 105);
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    char niz[8] = "Crtanje";
    BSP_LCD_DisplayStringAt(150, 160, (uint8_t*)niz, LEFT_MODE);
    char niz1[9] = "funkcije";
    BSP_LCD_DisplayStringAt(147, 190, (uint8_t*)niz1, LEFT_MODE);
}
void display_meni(){
    BSP_LCD_Clear(LCD_COLOR_BLUE);
    display_title();
    btn_settings();
    btn_interpolation();
    btn_functions();
}
void display_selected(double pot_read) {
    display_meni();
    int selection(pot_read * 3);
    switch(selection) {
        case 0: {
            BSP_LCD_SetTextColor(LCD_COLOR_RED);
            BSP_LCD_FillRect (9, 125, 110, 3);
            BSP_LCD_FillRect (9, 125, 3, 105);
            BSP_LCD_FillRect (9, 227, 110, 3);
            BSP_LCD_FillRect (116, 125, 3, 105);
            return;
        }
        case 1: {
            BSP_LCD_SetTextColor(LCD_COLOR_RED);
            BSP_LCD_FillRect (120, 125, 110, 3);
            BSP_LCD_FillRect (120, 125, 3, 105);
            BSP_LCD_FillRect (227, 125, 3, 105);
            BSP_LCD_FillRect (120, 227, 110, 3);
            return;
        }
        default: {
            BSP_LCD_SetTextColor(LCD_COLOR_RED);
            BSP_LCD_FillRect (200, 75, 30, 2);
            BSP_LCD_FillRect (200, 75, 2, 30);
            BSP_LCD_FillRect (200, 103, 30, 2);
            BSP_LCD_FillRect (228, 75, 2, 30);
        }
    }
}

// function selection

void display_functions_title(){
    BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
    BSP_LCD_SetFont(&Font24);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    char niz1[10] = "Odaberite";
    BSP_LCD_DisplayStringAt(0, 10, (uint8_t*)niz1, CENTER_MODE);
    char niz2[9] = "funkciju";
    BSP_LCD_DisplayStringAt(0, 40, (uint8_t*)niz2, CENTER_MODE);
    char niz3[11] = "za crtanje";
    BSP_LCD_DisplayStringAt(0, 70, (uint8_t*)niz3, CENTER_MODE);
}
void display_functions_selection(double pot_read){
    BSP_LCD_Clear(LCD_COLOR_BLUE);

    int i (pot_read * 12);
    
    display_functions_title();
    BSP_LCD_SetTextColor(LCD_COLOR_GRAY);
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    BSP_LCD_FillRect ( 10, 100, 220, 110);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    if (i == 0){
        //SINUS
        BSP_LCD_SetFont(&Font24);
        char niz[6] = "Sinus";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[15] = "a*sin(b*x+c)";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "1/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 1){
        //KOSINUS
        BSP_LCD_SetFont(&Font24);
        char niz[8] = "Cosinus";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char izraz[15] = "a*cos(b*x+c)";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)izraz, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "2/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 2){
        //SINC
        BSP_LCD_SetFont(&Font24);
        char niz[8] = "Sinc";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char izraz[15] = "a*sin(b*x+c)/x";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)izraz, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "3/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 3){
        //TESTERA
        BSP_LCD_SetFont(&Font24);
        char niz[8] = "Testera";
        BSP_LCD_DisplayStringAt(0, 140, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "4/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 4){
        //TROKUT
        BSP_LCD_SetFont(&Font24);
        char niz[8] = "Trokut";
        BSP_LCD_DisplayStringAt(0, 140, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "5/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 5){
        //CETVRTKA
        BSP_LCD_SetFont(&Font24);
        char niz[9] = "Cetvrtka";
        BSP_LCD_DisplayStringAt(0, 140, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "6/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 6){
        //EKSPONENCIJALNA
        BSP_LCD_SetFont(&Font20);
        char niz[17] = "Eksponencijalna";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[15] = "a*b^(x+c)";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "7/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 7){
        //LOGARITAMSKA
        BSP_LCD_SetFont(&Font24);
        char niz[17] = "Logaritamska";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[15] = "a*log(b*x)+c";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "8/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 8){
        //LINEARNA
        BSP_LCD_SetFont(&Font24);
        char niz[17] = "Linearna";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[15] = "a*x+b";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "9/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 9){
        //KVADRAtNA
        BSP_LCD_SetFont(&Font24);
        char niz[17] = "Kvadratna";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[15] = "a*x^2+b*x+c";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[6] = "10/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 10){
        //KUBNA
        BSP_LCD_SetFont(&Font24);
        char niz[17] = "Kubna";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[16] = "a*x^3+b*x^2+c*x";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[6] = "11/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else {
        //RECIPROCNA
        BSP_LCD_SetFont(&Font24);
        char niz[17] = "Reciprocna";
        BSP_LCD_DisplayStringAt(0, 130, (uint8_t*)niz, CENTER_MODE);
        BSP_LCD_SetFont(&Font20);
        char sinus[16] = "a/x+b";
        BSP_LCD_DisplayStringAt(0, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[6] = "12/12";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
}

// interpol selection

void display_interpolators_title(){
    BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
    BSP_LCD_SetFont(&Font24);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    char niz1[10] = "Odaberite";
    BSP_LCD_DisplayStringAt(0, 25, (uint8_t*)niz1, CENTER_MODE);
    char niz2[13] = "interpolator";
    BSP_LCD_DisplayStringAt(0, 55, (uint8_t*)niz2, CENTER_MODE);
}

void display_interpolators_selection(double pot_read){
    BSP_LCD_Clear(LCD_COLOR_BLUE);

    display_interpolators_title();

    int i (pot_read * 5);
    
    BSP_LCD_SetTextColor(LCD_COLOR_GRAY);
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    BSP_LCD_FillRect ( 10, 100, 220, 110);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    if (i == 0){
        BSP_LCD_SetFont(&Font24);
        char niz[9] = "Linearni";
        BSP_LCD_DisplayStringAt(0, 125, (uint8_t*)niz, CENTER_MODE);
        char sinus[15] = "interpolator";
        BSP_LCD_DisplayStringAt(0, 165, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "1/5";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 1){
       BSP_LCD_SetFont(&Font24);
        char niz[11] = "Polinomski";
        BSP_LCD_DisplayStringAt(0, 125, (uint8_t*)niz, CENTER_MODE);
        char sinus[15] = "interpolator";
        BSP_LCD_DisplayStringAt(0, 165, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "2/5";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 2){
        BSP_LCD_SetFont(&Font16);
        char niz[17] = "Trigonometrijski";
        BSP_LCD_DisplayStringAt(5, 130, (uint8_t*)niz, CENTER_MODE);
        char sinus[15] = "interpolator";
        BSP_LCD_DisplayStringAt(5, 170, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "3/5";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else if(i == 3){
        BSP_LCD_SetFont(&Font20);
        char niz[14] = "Baricentricni";
        BSP_LCD_DisplayStringAt(0, 125, (uint8_t*)niz, CENTER_MODE);
        char sinus[15] = "interpolator";
        BSP_LCD_DisplayStringAt(0, 165, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "4/5";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
    else {
        BSP_LCD_SetFont(&Font24);
        char niz[15] = "Spline";
        BSP_LCD_DisplayStringAt(0, 125, (uint8_t*)niz, CENTER_MODE);
        char sinus[15] = "interpolator";
        BSP_LCD_DisplayStringAt(0, 165, (uint8_t*)sinus, CENTER_MODE);
        BSP_LCD_SetFont(&Font16);
        BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
        BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
        char broj[5] = "5/5";
        BSP_LCD_DisplayStringAt(0, 220, (uint8_t*)broj, RIGHT_MODE);
    }
}

// settings meni

void display_settings_title(){
    BSP_LCD_SetBackColor(LCD_COLOR_BLUE);
    BSP_LCD_SetFont(&Font24);
    BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
    char niz1[10] = "Postavke";
    BSP_LCD_DisplayStringAt(0, 10, (uint8_t*)niz1, CENTER_MODE);
}
int pause_ms (0);
void display_settings_live(double potRead){
    char niz[5] = "Live";
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    BSP_LCD_DisplayStringAt(0, 50, (uint8_t*)niz, CENTER_MODE);
    char niz1[9] = "crtanje";
    BSP_LCD_DisplayStringAt(0, 75, (uint8_t*)niz1, CENTER_MODE);
    char niz2[11] = "Pauza[ms]";
    BSP_LCD_DisplayStringAt(0, 120, (uint8_t*)niz2, CENTER_MODE);
    pause_ms = 100 * potRead;

    char buf[5];
    sprintf(buf, "%d", pause_ms);
    BSP_LCD_DisplayStringAt(0, 145, (uint8_t*)buf, CENTER_MODE);
}
int colour (LCD_COLOR_WHITE);
void display_settings_colour(double potRead){
    char niz[5] = "Boja";
    BSP_LCD_SetBackColor(LCD_COLOR_GRAY);
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    BSP_LCD_DisplayStringAt(0, 50, (uint8_t*)niz, CENTER_MODE);
    char niz1[9] = "pozadine";
    BSP_LCD_DisplayStringAt(0, 75, (uint8_t*)niz1, CENTER_MODE);
    
    int code = potRead * 7;
    colour = LCD_COLOR_WHITE;
    if (code == 0) colour = LCD_COLOR_BLACK;
    else if (code == 1) colour = LCD_COLOR_BLUE;
    else if (code == 2) colour = LCD_COLOR_GRAY;
    else if (code == 3) colour = LCD_COLOR_RED;
    else if (code == 4) colour = LCD_COLOR_GREEN;
    else if (code == 5) colour = LCD_COLOR_YELLOW;

    
    BSP_LCD_SetTextColor(LCD_COLOR_BLACK);
    BSP_LCD_DrawRect(49, 109, 141, 81);
    //BSP_LCD_SetTextColor(LCD_COLOR_WHITE);
    BSP_LCD_SetTextColor(colour);
    BSP_LCD_FillRect ( 50, 110, 140, 79);
}

void display_settings(double potRead) {
    BSP_LCD_Clear(LCD_COLOR_BLUE);
    display_settings_title();
    BSP_LCD_SetTextColor(LCD_COLOR_GRAY);
    BSP_LCD_FillRect ( 30, 40, 180, 170);

    int setting(potRead * 2);
    if(setting == 0) display_settings_colour(B.read());
    else display_settings_live(B.read());
}


/////////////////////////////////////CONTROLER FUNCTIONS
double meni_selector(0);
double settings_slider(0);

void off() { 
    p.set_state(0);
    BSP_LCD_Clear(LCD_COLOR_BLACK); 
}

void meni() {
    p.set_state(1);
    p.clear_interpolate_dots();
    meni_selector = A.read();
    display_selected(meni_selector);
}

void f_selection() {
    p.set_state(2);
    meni_selector = A.read();
    display_functions_selection(meni_selector);
}
void draw_f() {
    p.set_state(3);
    p.get_value(X.read());
}
void dots_input() {
    p.set_state(4);
    p.inputInterpolateDots(false);
}
void i_selection() {
    p.set_state(5);
    meni_selector = A.read();
    display_interpolators_selection(meni_selector);
}
void draw_i() {
    p.set_state(6);
    p.get_value(X.read(), true);
}
void settings() {
    p.set_state(7);
    settings_slider = A.read();
    display_settings(settings_slider);
}

std::function<void()> event;


void onActivatePressed() {
    static bool on(true);
    on = !on;
    on? off() : meni();
    event = on? &off : &meni;
}

void onConfirmPressed() {
    unsigned short state (p.get_state());

    switch (state) {
        case 1: { //go to selected option(draw fun, interpolate or settings)
            int selected(meni_selector * 3);
            if(selected == 0) //interpole
                event = &dots_input;
            else if(selected == 1)  //draw fun
                event = &f_selection;
            else { //settings
                settings();
                event = &settings;
            }
            break;
        }
        case 2: {//draw selected function
            int selected(meni_selector * 12);
            switch(selected) {
                case 0: {
                    p.set_function(get_sin);
                    break;
                }
                case 1: {
                    p.set_function(get_cos);
                    break;
                }
                case 2: {
                    p.set_function(get_sinc);
                    break;
                }
                case 3: {
                    p.set_function(get_saw);
                    break;
                }
                case 4: {
                    p.set_function(get_triangle);
                    break;
                }
                case 5: {
                    p.set_function(get_squares);
                    break;
                }
                case 6: {
                    p.set_function(get_exponential);
                    break;
                }
                case 7: {
                    p.set_function(get_logarithmic);
                    break;
                }
                case 8: {
                    p.set_function(get_linearf);
                    break;
                }
                case 9: {
                    p.set_function(get_squaref);
                    break;
                }
                case 10: {
                    p.set_function(get_cubef);
                    break;
                }
                default: {
                    p.set_function(get_inversef);
                }
            }
            event = &draw_f;
            break;
        }
        case 3: {//doese nothing when function is drawn
            break;
        }
        case 4: {//interpole dots input
            i_selection();
            event = &i_selection;
            break;
        }
        case 5: {//interpolators selection
            int selected(meni_selector * 5);
            switch(selected) {
                case 0: {
                    lin = LinearInterpolator(p.get_interpolate_dots());
                    interpolator = &lin;
                    break;
                }
                case 1: {
                    pin = PiecewisePolynomialInterpolator(p.get_interpolate_dots(), 3);
                    interpolator = &pin;
                    break;
                }
                case 2: {
                    tin = TrigonometricInterpolator(p.get_interpolate_dots());
                    interpolator = &tin;
                    break;
                }
                case 3: {
                    bin = BarycentricInterpolator(p.get_interpolate_dots(), 3);
                    interpolator = &bin;
                    break;
                }
                default: {
                    spin = SplineInterpolator(p.get_interpolate_dots());
                    interpolator = &spin;
                    break;
                }
            }
            p.set_function(get_interpolate);
            draw_i();
            event = &draw_i;
           
        }
        case 6: {//draw interpole -> interolator selection
            break;
        }
        case 7: {
            
            int setting(settings_slider * 2);
            if(setting == 0) 
                p.set_background(colour);
            else 
                p.set_pause(pause_ms);
            break;
        }
        default: {
            //printf("ERROR!!!\n");
        }
    }
}

void onBackPressed() {
    unsigned short state (p.get_state());

    switch (state) {
        case 1: { //stays in main meni
            break;
        }
        case 2: {
            meni();
            event = &meni;
            break;
        }
        case 3: {
            f_selection();
            event = &f_selection;
            break;
        }
        case 4: {//interpole dots input
            meni();
            event = &meni;
            break;
        }
        case 5: {//interpolators selection
            dots_input();
            event = &dots_input;
            break;
        }
        case 6: {//draw interpole -> interolator selection
            i_selection();
            event = &i_selection;
            break;
        }
        case 7: {//settings -> main menu
            meni();
            event = &meni;
            break;
        }
        default: {
            //printf("ERROR!!!\n");
        }
    }
}

int main() {
    back_btn.rise(&onBackPressed);
    power_btn.rise(&onActivatePressed);
    confirm_btn.rise(&onConfirmPressed);
    p.screen_setup();

    event = &meni;
    for(;;) {
       event();
       wait_ms(100);        
    }
    return 0;
}
