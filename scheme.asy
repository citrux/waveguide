settings.outformat="pdf";
settings.inlineimage=true;
//settings.embed=true;
settings.toolbar=false;
viewportmargin=(2,2);

import patterns;
size(10cm);

add("hatch",hatch(1mm));

filldraw((-0.1,-0.1)--(5.1,-0.1)--(5.1,3.1)--(-0.1,3.1)--cycle,
pattern("hatch"));
filldraw((0,0)--(5,0)--(5,3)--(0,3)--cycle, white);
pen bold=linewidth(2*linewidth());
draw((1.5,0)--(1.5,3), bold);

draw((0,0)--(0,-1.2));
draw((1.5,0)--(1.5,-1.2));
draw((0,-1)--(1.5,-1), Arrows);
label("$c$", (0.75, -1), N);

draw((0,3)--(0,4), Arrow); label("$y$", (0,4), W);
draw((5,3)--(5,3.7));
draw((0,3.5)--(5,3.5), Arrows);
label("$a$", (2.5, 3.5), N);

draw((5,3)--(5.7,3));
draw((5,0)--(6,0), Arrow); label("$x$", (6.0), S);
draw((5.5,0)--(5.5,3), Arrows);
label("$b$", (5.5, 1.5), E);


label("$\varepsilon_2$", (0.75, 1.5));
label("$\varepsilon_1$", (3.25, 1.5));
