#set text(size: 8pt,font: "New Computer Modern",)
#set page(paper: "a4",flipped: true,margin: (x:1cm,y:1cm))
#set page(header: context [
#block(inset:5pt)[
  *University of Latvia*
  #h(1fr)
  #counter(page).display(
    "1/1",
    both: true,
  )
]
])

#set par(justify: true)
#set document(title: "LU ICPC kladīte ;)",author: ("Krišjānis Petručeņa","Matīss Kristiņš", "Valters Kalniņš", "Kristaps Štāls"))
#set heading(numbering: "1.")
// #show: columns.with(3, gutter: 2em)
#columns(3, gutter: 2em)[

#show heading.where(
  level: 1
): it => block(width: 100%)[
  #set align(center)
  #set text(12pt, weight: "regular")
  #it.body
  #v(1em)
]

#align(center)[#block(text(weight: 700, 1.75em, "LU ICPC kladīte ;)"))]
#outline(indent: 2em)

#include "cpplang.typ"

#colbreak()

#include "algebra.typ"

#colbreak()

#include "geometry.typ"

#colbreak()

#include "structures.typ"

#colbreak()

#include "graph.typ"

#colbreak()

#include "strings.typ"

#colbreak()

#include "dynamic.typ"

#colbreak()

]

#include "other.typ"