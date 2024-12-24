# Информация по содержанию программы    
## Модули прораммы
Astro.cpp - Основная функция - CalculateLatLonAlt. Используется для получения геодезических координат спутника. </br>
Math.cpp - Реализует различные математические методы.    
vector.h - Самописный вектор.  
time.cpp - Содержит методы для работы со временем.  
observer.cpp - Не используется.  
Utils.cpp - Содержит основные методы класса.  
Sgpsdp.cpp - Содержит основные астрономические расчёты.  
Sdpsgp.h - Описание Класса и некоторых структур.  
main.cpp - Содержит сам алгоритм. 
</br>
## Работа программы
Программе в одной директории необходим файл, содержащий TLE данные спутника - TLE.txt    
Программа принимает на вход геодезические координаты наземной станции связи (широта, долгота и высота над уровнем Земли). Ввод осуществляется через пробел или Enter. Есть ограничения широты(-90; 90) и долготы(-180; 180)        
Выводит: Номер спутника по NORAD - модуль скорости спутника - дату и время (в формате UTC/GMT+0) его появления в зоне видимости станции


## upgrade-1
1. В качестве параметров на входе иметь фильтр:
float minAzm, float maxAzm, float minElv, float maxElv, int timeMinObserveSec
где  Azm (0-360) Elv (0 - 360) параметры сектора Антенны, а timeObserveSec минимальное время нахождения спутника в указанной зоне в секундах. Т е, например:
56.14717 37.76275 100 // 120.0 200.0 10.0 80.0 50
2. отобразить толлько спутники, которые будут больше 50 секунд в секторе  аzm: (120-200) elv: (10-80), 
3. в выводе пометить "+" (на нас летит) "-" (от нас летит) - Comlplited
4. отобразить имя спутника (необязательно)

\left(z\cdot\cos\left(\beta\right)\ -\ \sin\left(\beta\right)\cdot\left(y\cdot\sin\left(\alpha\right)+x\cdot\cos\left(\alpha\right)\right)\right)^{2\ \ }+\left(y\cdot\cos\left(\alpha\right)-x\cdot\sin\left(\alpha\right)\right)^{2}-\left(z\sin(\beta)+\cos(\beta)(y\sin(\alpha)+x\cos(\alpha))-R\right)^{2}=0\left\{10>z\sin(\beta)+\cos(\beta)(y\sin(\alpha)+x\cos(\alpha))-R>\ 0\right\}
x = \left(z\cdot\cos\left(\beta\right)\ -\ \sin\left(\beta\right)\cdot\left(y\cdot\sin\left(\alpha\right)+x\cdot\cos\left(\alpha\right)\right)-60\right)
y = \left(y\cdot\cos\left(\alpha\right)-x\cdot\sin\left(\alpha\right)\right)
z = \left(z\sin(\beta)+\cos(\beta)(y\sin(\alpha)+x\cos(\alpha))-R\right)

f\left(x,y,z\right)=\left(z\cdot\cos\left(\beta\right)\ -\ \sin\left(\beta\right)\cdot\left(y\cdot\sin\left(\alpha\right)+x\cdot\cos\left(\alpha\right)\right)-60\right)

g\left(x,y,z\right)=\left(y\cdot\cos\left(\alpha\right)-x\cdot\sin\left(\alpha\right)\right)

h\left(x,y,z\right)=\left(z\sin(\beta)+\cos(\beta)(y\sin(\alpha)+x\cos(\alpha))-10\right)

\left(f\left(x,y,z\right)\cdot\cos\left(v\right)\ -\ \sin\left(v\right)\cdot\left(g\left(x,y,z\right)\cdot\sin\left(c\right)+h\left(x,y,z\right)\cdot\cos\left(c\right)\right)\right)^{2\ \ }+\left(g\left(x,y,z\right)\cdot\cos\left(c\right)-h\left(x,y,z\right)\cdot\sin\left(c\right)\right)^{2}-\left(f\left(x,y,z\right)\sin(v)+\cos(v)(g\left(x,y,z\right)\sin(c)+h\left(x,y,z\right)\cos(c))\right)^{2}=0\left\{10>f\left(x,y,z\right)\sin(v)+\cos(v)(f\left(x,y,z\right)\sin(c)+h\left(x,y,z\right)\cos(c))>\ 0\right\}