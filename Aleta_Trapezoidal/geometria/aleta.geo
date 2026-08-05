// =====================================================================
//  aleta.geo - Aleta trapezoidal com furos circulares internos
//  Geometria parametrizada para o subprojeto de IC (MFDM/MWLS)
//  UFES - Metodos Numericos em Transferencia de Calor
//
//  USO:
//    gmsh aleta.geo -2                          (parametros padrao)
//    gmsh -setnumber lc 0.10 aleta.geo -2        (refina a nuvem)
//    gmsh -setnumber r  0.30 aleta.geo -2        (muda o raio do furo)
// =====================================================================

SetFactory("OpenCASCADE");

// ---------------------------------------------------------------------
// 1. PARAMETROS  (DefineConstant permite sobrescrever por -setnumber)
// ---------------------------------------------------------------------
DefineConstant[ Wtopo = {6.0,  Name "Largura do topo"  } ];
DefineConstant[ Wbase = {4.0,  Name "Largura da base"  } ];
DefineConstant[ H     = {3.5,  Name "Altura total"     } ];
DefineConstant[ r     = {0.45, Name "Raio dos furos"   } ];
DefineConstant[ xf    = {1.5,  Name "Offset x do furo" } ];
DefineConstant[ yf    = {1.75, Name "Altura y do furo" } ];
DefineConstant[ lc    = {0.15, Name "Tamanho de malha global" } ];
DefineConstant[ fref  = {0.40, Name "Fator de refino na borda do furo" } ];

tol = 1e-3;   // tolerancia para selecao por bounding box

// ---------------------------------------------------------------------
// 2. CONTORNO EXTERNO (trapezio)
//    Base menor embaixo (Wbase), base maior em cima (Wtopo)
// ---------------------------------------------------------------------
Point(1) = {-Wbase/2, 0, 0, lc};
Point(2) = { Wbase/2, 0, 0, lc};
Point(3) = { Wtopo/2, H, 0, lc};
Point(4) = {-Wtopo/2, H, 0, lc};

Line(1) = {1, 2};   // base   (inferior)
Line(2) = {2, 3};   // lado direito  (inclinado)
Line(3) = {3, 4};   // topo   (superior)
Line(4) = {4, 1};   // lado esquerdo (inclinado)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ---------------------------------------------------------------------
// 3. FUROS (discos completos - so o OpenCASCADE faz isso em 1 comando)
// ---------------------------------------------------------------------
Disk(2) = {-xf, yf, 0, r};
Disk(3) = { xf, yf, 0, r};

// ---------------------------------------------------------------------
// 4. SUBTRACAO BOOLEANA: aleta = trapezio - furos
//    Delete remove as entidades originais apos a operacao
// ---------------------------------------------------------------------
BooleanDifference(4) = { Surface{1}; Delete; }{ Surface{2, 3}; Delete; };

// ---------------------------------------------------------------------
// 5. IDENTIFICACAO DAS FRONTEIRAS
//    Apos a booleana as tags sao renumeradas -> selecionar por posicao,
//    nunca por numero fixo. "In BoundingBox" exige a curva INTEIRA dentro.
// ---------------------------------------------------------------------
cBase[]  = Curve In BoundingBox { -Wbase/2-tol, -tol,   -tol,
                                   Wbase/2+tol,  tol,    tol };

cTopo[]  = Curve In BoundingBox { -Wtopo/2-tol,  H-tol, -tol,
                                   Wtopo/2+tol,  H+tol,  tol };

cLadoD[] = Curve In BoundingBox {  Wbase/2-tol,  -tol,  -tol,
                                   Wtopo/2+tol,  H+tol,  tol };

cLadoE[] = Curve In BoundingBox { -Wtopo/2-tol,  -tol,  -tol,
                                  -Wbase/2+tol,  H+tol,  tol };

cFuroE[] = Curve In BoundingBox { -xf-r-tol, yf-r-tol, -tol,
                                  -xf+r+tol, yf+r+tol,  tol };

cFuroD[] = Curve In BoundingBox {  xf-r-tol, yf-r-tol, -tol,
                                   xf+r+tol, yf+r+tol,  tol };

// ---------------------------------------------------------------------
// 6. GRUPOS FISICOS
//    O NUMERO da tag vira o codigo de BC no arquivo exportado:
//      1 = topo  | 2 = base | 3 = lados | 4 = furos | 10 = dominio
// ---------------------------------------------------------------------
Physical Curve("topo",  1) = { cTopo[]  };
Physical Curve("base",  2) = { cBase[]  };
Physical Curve("lados", 3) = { cLadoD[], cLadoE[] };
Physical Curve("furos", 4) = { cFuroE[], cFuroD[] };
Physical Surface("aleta", 10) = {4};

// ---------------------------------------------------------------------
// 7. CAMPO DE REFINO NA BORDA DOS FUROS
//    Essencial para o item (iii): a concentracao de fluxo acontece na
//    borda do furo; sem refino local o pico de |q| fica subestimado.
// ---------------------------------------------------------------------
Field[1] = Distance;
Field[1].CurvesList = { cFuroE[], cFuroD[] };
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc * fref;   // malha fina colada no furo
Field[2].SizeMax  = lc;          // malha global longe do furo
Field[2].DistMin  = r * 0.5;     // ate 0.5r do furo -> SizeMin
Field[2].DistMax  = r * 3.0;     // alem de 3r      -> SizeMax

Background Field = 2;

// ---------------------------------------------------------------------
// 8. OPCOES DE MALHA
// ---------------------------------------------------------------------
Mesh.Algorithm = 6;              // Frontal-Delaunay (triangulos regulares)
Mesh.MshFileVersion = 2.2;       // ASCII 2.2 = facil de ler no Octave
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
