# nurbs-tessellation
Creating a nurbs surface and converting it to a mesh of triangles (tessellation).

Esse projeto foi feito em C++ e tem como biblioteca básica SFML no Linux (Ubuntu 22.04). O projeto implementado é "1.Superfície B-Spline Racional (NURBS) (Diego)" descrito em MG-1-2024-Projetos_.pdf.

Veja a documentação de como instalar pacotes de SFML em "https://www.sfml-dev.org/tutorials/2.6/start-linux.php".

*Para compilar código objeto:

g++ nurbs_tessellation02.cpp -o tessel -lsfml-graphics -lsfml-window -lsfml-system

*Para executar:

./tessel < in_nurbs.data > out_nurbs.data
