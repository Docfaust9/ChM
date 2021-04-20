#include <iostream>
#include<vector>
#include<cmath>
#include<iomanip>
using namespace std;
class Number
{
private:
    double count(double x);
    double phi(double x);
    double derivative(double x);
    double derivative2(double x);
    void find_borders(int x);
    double borderline;
public:
    int counter,root_counter;
    double a, b, epsi;
    vector<double> xpow, xk;
    void dihotom();
    void simple_iter();
    void sek();
    void newton();
    void hord();
    void init();
};
void Number::find_borders(int x)
{
    a = x;
    for (int i = x+1; i <= 1000; i++)
    {
        b = i;
        if (count(a) * count(b) > 0)
        {
            a = i;
        }
        else
            if (count(a) * count(b) < 0)
            {
                break;
            }
     
    }
    borderline = b;
}
double Number::phi(double x)
{
    double result;
    double buf = 0;
    for (int i = 1; i < counter; i++)
    {
        buf += -xk[i] * pow(x, xpow[i]);
    }
    result = pow(buf, 1 / xpow[0]) / xk[0];
    return result;
}
double Number::derivative(double x)
{
    double result = 0;
    for (int i = 0; i < counter; i++)
    {
        result += xk[i]*xpow[i] * pow(x, (xpow[i]-1));
    }
    return result;

}
double Number::derivative2(double x)
{
    double result = 0;
    for (int i = 0; i < counter; i++)
    {
        result += xk[i] * xpow[i] *(xpow[i]-1)* pow(x, (xpow[i]-2));
    }
    return result;

}
double Number::count(double x) 
{
    double result=0;
    for (int i = 0; i < counter; i++)
    {
        result += xk[i] * pow(x, xpow[i]);
    }
    return result;
}
void Number::init()
{
    cout << "Введите погрешность(epsilon): ";
    cin >> epsi;
    cout << "Введите количество слагаемых: ";
    cin >> counter;
    xk.resize(counter);
    xpow.resize(counter);
    cout << "Введите коэффиценты при каждом слагаемоем: "<<endl;
    for (int i = 0; i < counter; i++)
    {
        cout << i + 1 << "-е слагаемое: ";
        cin >> xk[i];
    }
    cout << "Введите степень х при каждом из слагаемых: "<<endl;
    for (int i = 0; i < counter; i++)
    {
        cout << i + 1 << "-е слагаемое: ";
        cin >> xpow[i];
    }
    cout << endl<<"Вы ввели: ";
    for (int i = 0; i <counter; i++)
    {
        if ((xk[i] >= 0) && (i != 0))
            cout << "+";
        if ((xpow[i] != 0)&&(xpow[i]!=1))
        {
            cout << xk[i]  << "*x^" << xpow[i];
        }
        else
            if(xpow[i]==0)
            {
                cout << xk[i];
            }
            else
                if (xpow[i] == 1)
                {
                    cout << xk[i]<<"*x";
                }
    }
    cout << "=0"<<endl<<endl;
    cout << "Укажите количество действительных корней: ";
    cin >> root_counter;
    cout << endl;
}
void Number::dihotom()
{
    borderline = -1000;
    cout << "Методом половинного деления: " << endl;
    for (int j = 0; j < root_counter; j++)
    {
        find_borders(borderline);
        double x;
        while (abs(a - b) > (2 * epsi))
        {
            x = (a + b) / 2;
            if (count(a) * count(x) < 0)
                b = x;
            else
                a = x;
        }
        x = (a + b) / 2;
        cout << "x" << j + 1 << "=" << setprecision(5)<< x<<endl;
    }
}
void Number::simple_iter()
{
    borderline = -1000;
    cout << "Методом простых итераций: " << endl;
    for (int j = 0; j < root_counter; j++)
    {
        find_borders(borderline);
        double x0 = (a + b) / 2;
        double x = phi(x0);
        while (abs(x - x0) > epsi)
        {
            x0 = x;
            x = phi(x0);
        }
        cout << "x" << j + 1 << "=" << setprecision(5) << x << endl;
    }
}
void Number::sek()
{
    borderline = -1000;
    cout << "Методом секущих: " << endl;
    for (int j = 0; j < root_counter; j++)
    {
        find_borders(borderline);
        double x0 = a;
        double x1 = b;

        while (abs(x1 - x0) > epsi)
        {
            double tmp = x1;
            x1 = x1 - ((x1 - x0) * count(x1) / (count(x1) - count(x0)));
            x0 = tmp;
        }
        cout << "x" << j + 1 << "=" << setprecision(5) << x1 << endl;
    }

}
void Number::newton()
{
    borderline = -1000;
    cout << "Методом Ньютона(касательных): " << endl;
    for (int j = 0; j < root_counter; j++)
    {
        find_borders(borderline);
        double delta, x;
        if ((count(a) * derivative2(a)) >= 0)
            x = a;
        else
            x = b;
        delta = count(x) / derivative(x);
        while (abs(delta) > epsi)
        {
            delta = count(x) / derivative(x);
            x = x - delta;
        }
        cout << "x" << j + 1 << "=" << setprecision(5) << x << endl;
    }
}
void Number::hord()
{
    borderline = -1000;
    cout << "Методом хорд(ложного положения): " << endl;
    for (int j = 0; j < root_counter; j++)
    {
        find_borders(borderline);
        double x0 = a;
        double x1 = b;

        while (abs(x1 - x0) > epsi)
        {
            double tmp = x1;
            x1 = x0 - ((x1 - x0) * count(x0) / (count(x1) - count(x0)));
            x0 = tmp;
        }
        cout << "x" << j + 1 << "=" << setprecision(5) << x1 << endl;
    }

}
int main()
{
    setlocale(LC_ALL, "Russian");
    Number obj;
    obj.init();
    cout << fixed;
    int choose;
    cout << "Выберите численный метод для поиска корня(корней)" << endl << "1. Метод половинного деления\n2. Метод простых итераций\n3. Метод секущих\n4. Метод Ньютона(касательных)\n5. Метод хорд(ложного положения)\n0. Выход"<<endl;
    cin >> choose;
    while (choose != 0)
    {
        switch (choose)
        {
        case 1:
            obj.dihotom();
            break;
        case 2:
            obj.simple_iter();
            break;
        case 3:
            obj.sek();
            break;
        case 4:
            obj.newton();
            break;
        case 5:
            obj.hord();
            break;
        default: cout<<"Неверное значение!"<<endl;
        }
        cin >> choose;
    }
    return 0;
}


