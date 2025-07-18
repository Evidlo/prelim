// This function gets your whole document as its `body` and formats
// it as an article in the style of the IEEE.
#let ieee(
    // The paper's title.
    title: "Paper Title",

    // An array of authors. For each author you can specify a name,
    // department, organization, location, and email. Everything but
    // but the name is optional.
    authors: (),

    // The paper's abstract. Can be omitted if you don't have one.
    abstract: none,

    // A list of index terms to display after the abstract.
    index-terms: (),

    // The article's paper size. Also affects the margins.
    paper-size: "us-letter",

    // The path to a bibliography file if you want to cite some external
    // works.
    bibliography-file: none,

    // The paper's content.
    body
) = {
    // Set document metdata.
    set document(title: title, author: authors.map(author => author.name))

    // Set the body font.
    set text(font: "New Computer Modern", size: 12pt)

    // Configure the page.
    set page(
        // paper: paper-size,
        // FIXME - for Lara's notes
        width: 12.5in, height: 11in,
        // The margins depend on the paper size.
        margin: (
            x: 1.45in,
            top: 5%,
            bottom: 5%,
            // FIXME - for Lara's notes
            right: 5.0119in
        ),
        numbering: "1"
    )

    // Configure equation numbering and spacing.
    show math.equation: set block(spacing: 0.65em)
    set math.equation(numbering: "(1)", block: true)
    // only show number on equations with labels
    show math.equation: it => {
        if it.block and not it.has("label") [
            // dont let unlabeled equations increment number
            #counter(math.equation).update(n => calc.max(n - 1, 0))
            #math.equation(it.body, block: true, numbering: none)#label("")

        ] else {
            it
        }
    }

    // Give figures some spacing
    // FIXME: this is quite ugly
    // let figure_spacing = 1em // Additional spacing between figures and the text
    // show figure: it => {
    //     if it.placement == none {
    //         block(it, inset: (y: figure_spacing))
    //     } else if it.placement == top {
    //         place(
    //             it.placement,
    //             float: true,
    //             block(width: 100%, inset: (bottom: figure_spacing), align(center, it))
    //         )
    //     } else if it.placement == bottom {
    //         place(
    //             it.placement,
    //             float: true,
    //             block(width: 100%, inset: (top: figure_spacing), align(center, it))
    //         )
    //     }
    // }


    // Configure lists.
    set enum(indent: 10pt, body-indent: 9pt)
    set list(indent: 10pt, body-indent: 9pt)

    // Configure headings.
    set heading(numbering: "1.1.1.")
    show heading.where(level: 2): set text(14pt)
    // show heading: it => locate(loc => {
    // Find out the final number of the heading counter.
    //   let levels = counter(heading).at(loc)
    //   let deepest = if levels != () {
    //     counter(heading)
    //   } else {
    //     1
    //   }

    //   set text(10pt, weight: 400)
    //   if it.level == 1 [
    //     // First-level headings are centered smallcaps.
    //     // We don't want to number of the acknowledgment section.
    //     #let is-ack = it.body in ([Acknowledgment], [Acknowledgement])
    //     #set align(center)
    //     #set text(if is-ack { 10pt } else { 12pt })
    //     #show: smallcaps
    //     #v(20pt, weak: true)
    //     #if it.numbering != none and not is-ack {
    //       numbering("1.", deepest)
    //       h(7pt, weak: true)
    //     }
    //     #it.body
    //     #v(13.75pt, weak: true)
    //   ] else if it.level == 2 [
    //     // Second-level headings are run-ins.
    //     #set par(first-line-indent: 0pt)
    //     #set text(style: "italic")
    //     #v(10pt, weak: true)
    //     #if it.numbering != none {
    //       numbering("1.1.", deepest)
    //       h(7pt, weak: true)
    //     }
    //     #it.body
    //     #v(10pt, weak: true)
    //   ] else [
    //     // Third level headings are run-ins too, but different.
    //     #if it.level == 3 {
    //       numbering("1.1.1.", deepest)
    //       [ ]
    //     }
    //     _#(it.body):_
    //   ]
    // })

    // Display the paper's title.
    v(3pt, weak: true)
    align(center, text(26pt, title))
    v(8.35mm, weak: true)

    // Display the authors list.
    for i in range(calc.ceil(authors.len() / 3)) {
        let end = calc.min((i + 1) * 3, authors.len())
        let is-last = authors.len() == end
        let slice = authors.slice(i * 3, end)
        grid(
            columns: slice.len() * (1fr,),
            gutter: 12pt,
            ..slice.map(author => align(center, {
                text(12pt, author.name)
                if "department" in author [
                    \ #emph(author.department)
                ]
                if "organization" in author [
                    \ #emph(author.organization)
                ]
                if "location" in author [
                    \ #author.location
                ]
                if "email" in author [
                    \ #link("mailto:" + author.email)
                ]
            }))
        )

        if not is-last {
            v(16pt, weak: true)
        }
    }
    align(center, text(12pt, emph("Electrical and Computer Engineering\nUniversity of Illinois Urbana-Champaign")))
    v(40pt, weak: true)

    // Start two column mode and configure paragraph properties.
    show: columns.with(1, gutter: 0pt)
    set par(justify: true, first-line-indent: 1em, spacing:.8em, leading:.8em)

    // Display abstract and index terms.
    if abstract != none [
        #set text(weight: 700)
        #h(1em) _Abstract_---#abstract

        #if index-terms != () [
            #h(1em)_Index terms_---#index-terms.join(", ")
        ]
        #v(2pt)
    ]

    // Display the paper's contents.
    body

    // Display bibliography.
    if bibliography-file != none {
        show bibliography: set text(8pt)
        bibliography(bibliography-file, title: text(10pt)[References], style: "ieee")
    }
}
