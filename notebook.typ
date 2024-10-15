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

#colbreak()

#include "1-cpplang.typ"

#colbreak()

#include "2-algebra.typ"

#colbreak()

#include "3-geometry.typ"

#colbreak()

#include "4-structures.typ"

#colbreak()

#include "5-graph.typ"

#colbreak()

#include "6-strings.typ"

#colbreak()

#include "7-dynamic.typ"

#colbreak()

]

#include "8-other.typ"