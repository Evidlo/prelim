// arrow above vector variable
#let vc(symbol) = {
    $limits(symbol)^harpoon.rt$
}


// shaded box around content
#let rbox(content) = {
    rect(
        inset: .5em,
        // radius: .5em,
        fill: luma(245),
        content
    )
}

// vertical text - e.g. table header row
#let vert(content) = {
    rotate(content, -90deg, reflow: true)
}

// expectation
#let ex(content) = {
    $bb(E)[content]$
}
// variance
#let var(content) = {
    $"Var"[content]$
}

// inner product with brackets
#let innerproduct(x, y) = $lr(angle.l #x, #y angle.r)$

// subfigure with labelling
#let subfigure(content, kind, caption) = {
    set figure.caption(separator: " ")
    figure(content, caption: caption, numbering: "(a)", supplement: "", kind: kind)
}

// gray text - e.g. comments
#let gt(content) = {
    set text(fill: gray)
    content
}

// red text - e.g. draft text
#let rt(content) = {
    set text(fill: red)
    content
}

// gray text background
#let ga(content) = {
    highlight(fill: gray, extent: 3pt, radius: 8pt, content)
}

// showy box definition
#import "@preview/showybox:2.0.4": showybox
#let sb(title, content) = {
    showybox(
        title-style: (weight: "bold", boxed-style: (:)),
        title: title,
        content
    )
}

#let s(title, content) = {
    showybox(
        title-style: (boxed-style: (:)),
        title: title,
        content
    )
}