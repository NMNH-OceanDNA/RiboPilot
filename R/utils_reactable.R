# reactable rendering helpers (mirrors MitoPilot's utils_reactable.R).

#' Highlight selected row(s)
#'
#' Returns nothing for unselected rows so `striped = TRUE` still bands the
#' table; forcing #FFFFFF here silently cancelled it.
#' @noRd
rt_highlight_row <- function() {
  htmlwidgets::JS("
    function(rowInfo) {
      if (typeof rowInfo === 'undefined') return
      if (!rowInfo.selected) return null
      return { background: '#D3BEC2' }
    }")
}

#' Dynamic icon cell: map a cell value to a FontAwesome icon
#' @param icons named character vector (value -> "fa fa-..." class)
#' @noRd
rt_dynamicIcon <- function(icons = NULL) {
  if (length(icons) == 0) {
    return(htmlwidgets::JS("function(cellInfo) {return cellInfo.value}"))
  }
  key <- paste(paste0(names(icons), ": '", icons, "'"), collapse = ", ")
  sprintf(
    "function(cellInfo){
       var { value } = cellInfo;
       var icons = { %s };
       var icon = icons[value] || '';
       return `<i class='${icon}' style='padding-left: 0.2em;'></i>`
     }", key
  ) |> htmlwidgets::JS()
}

#' Column header with a "?" that explains the column on hover
#'
#' Use with `colDef(html = TRUE)`. The explanation rides on the native `title`
#' attribute, so it needs no JS and works wherever the widget is rendered.
#' @param name Header label. @param tip One-sentence explanation.
#' @noRd
rt_help_header <- function(name, tip) {
  as.character(shiny::tags$span(
    title = tip, style = "cursor: help;",
    name,
    shiny::tags$span(style = "color:#999; font-size:0.85em; margin-left:3px;", "?")))
}

#' Hover tooltip for truncated text (use with wrap = FALSE)
#' @noRd
rt_longtext <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var text = cellInfo.value ? cellInfo.value : ''
      return `<abbr style='cursor: info; text-decoration: none;' title='${text}'>${text}</abbr>`
    }")
}

#' Render a yes/no text column as a colored badge
#' @param invert if TRUE, 'yes' is orange (bad) and 'no' is green
#' @param hide_no if TRUE, render nothing for 'no'/NA
#' @noRd
rt_bool_badge <- function(invert = FALSE, hide_no = FALSE) {
  yes_bg <- if (invert) "#fde8d0" else "#d4edda"
  yes_fg <- if (invert) "#7d4a1e" else "#2d6a4f"
  no_bg <- if (invert) "#d4edda" else "#fde8d0"
  no_fg <- if (invert) "#2d6a4f" else "#7d4a1e"
  sprintf(
    "function(cellInfo) {
      var val = cellInfo.value ? cellInfo.value : 'no'
      if (val !== 'yes' && %s) return ''
      var bg = val === 'yes' ? '%s' : '%s'
      var fg = val === 'yes' ? '%s' : '%s'
      return '<span style=\"background:' + bg + '; color:' + fg + '; border-radius:3px; padding:1px 6px; font-size:0.85em;\">' + val + '</span>'
    }", tolower(as.character(hide_no)), yes_bg, no_bg, yes_fg, no_fg
  ) |> htmlwidgets::JS()
}

#' Render a curation flag count: blank/green when 0, orange badge when > 0
#' @noRd
rt_flag_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var n = cellInfo.value
      if (n === null || n === undefined || n === '' ) return ''
      n = Number(n)
      if (n === 0) return '<span style=\"color:#2d6a4f;\">&#10003;</span>'
      return '<span style=\"background:#fde8d0; color:#7d4a1e; border-radius:3px; padding:1px 6px; font-size:0.85em;\">' + n + '</span>'
    }")
}

#' Render a kept-scaffold count: red when nothing survives triage, else plain
#'
#' Zero kept units means the sample has no assembly to annotate or export, which
#' is otherwise indistinguishable from a normal low count.
#' @noRd
rt_kept_count <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var n = cellInfo.value
      if (n === null || n === undefined || n === '') return ''
      n = Number(n)
      if (n === 0) return '<span style=\"background:#f8d7da; color:#842029; border-radius:3px; padding:1px 6px; font-size:0.85em;\">0</span>'
      return String(n)
    }")
}

#' Render an NCBI accession as a clickable link
#' @noRd
rt_ncbi_link <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var text = cellInfo.value ? cellInfo.value : ''
      if (!text || text === 'NA') return text
      var url = 'https://www.ncbi.nlm.nih.gov/nuccore/' + text
      return `<a href='${url}' target='_blank' rel='noopener noreferrer'>${text}</a>`
    }")
}

#' Clickable cell that fires a Shiny input with the row id (index+1)
#' @param InputId shiny input id to set on click
#' @noRd
rt_link <- function(InputId) {
  sprintf(
    "function(cellInfo) {
      var clickid = '%s';
      var sampid = cellInfo.index + 1;
      return `<a href='#' id=${sampid} class='grow' ` +
        `onclick='event.stopPropagation(); Shiny.setInputValue(&#39;${clickid}&#39;, this.id, {priority: &#39;event&#39;})'>` +
        (cellInfo.value ? cellInfo.value : '') +
        `</a>`;
    }", InputId
  ) |> htmlwidgets::JS()
}

#' Tri-state keep / drop control, emitting "<uid>|keep" or "<uid>|drop"
#'
#' The third state is the ABSENCE of a decision, not a third button: a unit with
#' no verdict shows both outlines unfilled and an amber cast. Clicking the state
#' a row is already in clears it back to undecided, so a decision is reversible
#' without a separate control.
#' @param InputId shiny input id to set on click.
#' @noRd
rt_triage_bttn <- function(InputId) {
  sprintf(
    "function(cellInfo) {
      var clickid = '%s';
      var uid = cellInfo.row['uid'];
      var v = cellInfo.value;
      var reason = cellInfo.row['triage_reason'] ? cellInfo.row['triage_reason'] : '';
      var src = cellInfo.row['triage_source'];
      var keep = (v === 'keep'), drop = (v === 'drop');
      var tip = reason ? reason : (keep ? 'kept' : (drop ? 'dropped' : 'not yet decided'));
      if (src) tip += ' (' + src + ')';
      // Emits the full verdict word, not the button letter: the server side
      // takes it as a verdict, and a bare 'k' only resolved to 'keep' by
      // match.arg()'s partial matching while surfacing as \"marked k\" to the user.
      function btn(label, verdict, on, colOn) {
        var st = on
          ? 'background:' + colOn + '; color:#fff; border:1px solid ' + colOn + ';'
          : 'background:transparent; color:#999; border:1px solid ' + (v ? '#ccc' : '#e0a800') + ';';
        return '<a href=\"#\" title=\"mark ' + verdict + '\" ' +
          'style=\"display:inline-block; width:20px; text-align:center; ' +
          'font-size:0.75em; line-height:16px; border-radius:3px; margin-right:2px; ' +
          'text-decoration:none; ' + st + '\" ' +
          'onclick=\"event.stopPropagation(); Shiny.setInputValue(&#39;' + clickid +
          '&#39;, &#39;' + uid + '|' + verdict + '&#39;, {priority: &#39;event&#39;})\">' +
          label + '</a>';
      }
      return '<span title=\"' + tip.replace(/\"/g, '') + '\">' +
        btn('K', 'keep', keep, '#2d6a4f') + btn('D', 'drop', drop, '#c0392b') + '</span>';
    }", InputId) |> htmlwidgets::JS()
}

#' Render the deepest shared taxonomic rank between a hit and the sample
#'
#' Deliberately never asserts "host". The evidence supports "shares X at rank Y",
#' and a green host badge on a co-extracted copepod would be worse than showing
#' nothing. `no reference` renders NEUTRAL GREY, not amber: for an offline
#' project or a sparse clade it is the normal state, and painting it as
#' needs-attention makes every sample in the tool's core use case look broken.
#' @noRd
rt_shares_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === null || v === undefined || v === '' || v === 'NA')
        return '<span style=\"color:#999;\" title=\"no reference resolved for this scaffold\">no reference</span>'
      var close = {species:1, genus:1, family:1}
      var mid   = {order:1, class:1}
      var bg, fg
      if (close[v])      { bg = '#d4edda'; fg = '#2d6a4f' }
      else if (mid[v])   { bg = '#eef7e6'; fg = '#40693f' }
      else if (v === 'none') { bg = '#f8d7da'; fg = '#842029' }
      else               { bg = '#fde8d0'; fg = '#7d4a1e' }
      var t = cellInfo.row['screen_species'] ? cellInfo.row['screen_species'] : ''
      return '<span title=\"top GenBank hit: ' + t + '\" style=\"background:' + bg +
        '; color:' + fg + '; border-radius:3px; padding:1px 6px; font-size:0.85em;\">' + v + '</span>'
    }")
}

#' Which BLAST named the species in the neighbouring hit columns
#'
#' 'core_nt' is the remote GenBank search; 'local_rrna' is a local search
#' against a configured rDNA reference panel, a real reference search just run
#' locally instead of remotely; 'seeds' is the local search against the
#' GetOrganelle seed file, which is a hint rather than an identification. Rows
#' written before the source was recorded show as unknown, not as any of these.
#' @noRd
rt_blast_source_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === 'core_nt')
        return '<span title=\"remote megablast against NCBI core_nt\" style=\"background:#d4edda;' +
          ' color:#2d6a4f; border-radius:3px; padding:1px 6px; font-size:0.85em;\">GenBank</span>'
      if (v === 'local_rrna')
        return '<span title=\"local BLAST against the configured rDNA reference panel\"' +
          ' style=\"background:#d4edda; color:#2d6a4f; border-radius:3px; padding:1px 6px;' +
          ' font-size:0.85em;\">local panel</span>'
      if (v === 'seeds')
        return '<span title=\"local BLAST against the GetOrganelle seed file: the closest bait' +
          ' sequence, NOT a GenBank identification\" style=\"background:#fde8d0; color:#7d4a1e;' +
          ' border-radius:3px; padding:1px 6px; font-size:0.85em;\">seed file</span>'
      if (v === null || v === undefined || v === '' || v === 'NA') {
        var sp = cellInfo.row['blast_species']
        if (sp === null || sp === undefined || sp === '') return ''
        return '<span style=\"color:#999;\" title=\"assembled before the search was recorded\">unknown</span>'
      }
      return v
    }")
}

#' Flag a scaffold whose top hit disagrees with its own siblings
#' @noRd
rt_odd_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      if (cellInfo.value !== true && cellInfo.value !== 'TRUE') return ''
      var sib = cellInfo.row['sib_taxon'] ? cellInfo.row['sib_taxon'] : 'the others'
      var rk  = cellInfo.row['sib_rank'] ? cellInfo.row['sib_rank'] : ''
      return '<span title=\"this scaffold&#39;s top hit differs from the rest of the sample (' +
        sib + ' at ' + rk + ')\" style=\"background:#fde8d0; color:#7d4a1e; border-radius:3px;' +
        ' padding:1px 6px; font-size:0.85em;\">odd</span>'
    }")
}

#' Render the contamination-screen verdict: is the scaffold on this sample's
#' own lineage.
#'
#' `no_reference` renders NEUTRAL GREY like [rt_shares_badge()]'s "no
#' reference": for an offline project or a sparse clade it is the normal
#' state, not a problem. A missing screen row (never run) reads the same way,
#' since the app should not tell the two apart.
#' @noRd
rt_screen_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === null || v === undefined || v === '' || v === 'NA' || v === 'no_reference')
        return '<span style=\"color:#999;\">no_reference</span>'
      var bg, fg
      if (v === 'on_target')       { bg = '#d4edda'; fg = '#2d6a4f' }
      else if (v === 'off_target') { bg = '#f8d7da'; fg = '#842029' }
      else                          { bg = '#fde8d0'; fg = '#7d4a1e' }
      return '<span style=\"background:' + bg + '; color:' + fg +
        '; border-radius:3px; padding:1px 6px; font-size:0.85em;\">' + v + '</span>'
    }")
}

#' Render the chimera check verdict: clean, suspect, or likely, with the
#' reason as a hover tooltip.
#' @noRd
rt_chimera_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (v === null || v === undefined || v === '' || v === 'NA') return '';
      var bg, fg;
      if (v === 'clean')        { bg = '#d4edda'; fg = '#2d6a4f'; }
      else if (v === 'likely')  { bg = '#f8d7da'; fg = '#842029'; }
      else                      { bg = '#fde8d0'; fg = '#7d4a1e'; }
      var tip = cellInfo.row.chimera_reason || '';
      return '<span title=\"' + String(tip).replace(/\"/g, '&quot;') + '\" style=\"background:' + bg +
        '; color:' + fg + '; padding:1px 6px; border-radius:8px; font-size:0.85em\">' + v + '</span>';
    }")
}

#' Render the intragenomic variant class: is this scaffold one of several rDNA
#' copies, or the same molecule assembled twice.
#'
#' NEITHER class is a defect, so neither may borrow [rt_screen_badge()]'s
#' off-target red. An isomer renders NEUTRAL, because it is telling the curator a
#' decision does NOT matter. A paralog renders as INFORMATION: it is real
#' biology, documented in Octocorallia, Scyphozoa and other groups, and colouring
#' it like a problem is what this sub-project exists to stop. A `singleton` and a
#' unit with no finding render empty rather than as a badge, so a project that
#' has never been screened does not sprout a column of placeholders.
#'
#' The measured identity to the group representative goes in the title attribute,
#' so the number is reachable without widening the table.
#' @noRd
rt_variant_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === null || v === undefined || v === '' || v === 'NA' || v === 'singleton')
        return ''
      var row = cellInfo.row || {}
      var grp = row['variant_group']
      var pid = row['variant_identity']
      var tip = v === 'isomer'
        ? 'same molecule assembled two ways; differs only outside the conserved anchors, pick either'
        : 'this sample carries several rDNA copies; which you lock determines what GenBank receives'
      if (pid !== null && pid !== undefined && pid !== '')
        tip += ' (anchors ' + pid + ' percent identical to the group representative)'
      var bg, fg
      if (v === 'isomer') { bg = '#eeeeee'; fg = '#555555' }
      else                { bg = '#dbe9f6'; fg = '#1d4e79' }
      return '<span title=\"' + tip + '\" style=\"background:' + bg + '; color:' + fg +
        '; border-radius:3px; padding:1px 6px; font-size:0.85em;\">' + v +
        (grp ? ' ' + grp : '') + '</span>'
    }")
}

#' Inline depth bar relative to the deepest scaffold of the same sample
#'
#' The value is the MEDIAN per-position depth where one is stored (see
#' R/coverage.R for why a mean misleads on an rDNA scaffold); the tooltip names
#' which statistic it is, so a pre-v10 row falling back to the mean is not read
#' as a median.
#'
#' The bar has no background track. A grey remainder said only "this is a
#' fraction of something", which the fill length and the number beside it
#' already say. The 48px slot it used to shade is kept, transparent: the width
#' is what puts every row's number in the same place, and a bar that sets its
#' own width drags the number left and right by up to 46px down the column.
#' @param max_col Column holding the sample maximum.
#' @noRd
rt_depth_bar <- function(max_col = "depth_max") {
  sprintf(
    "function(cellInfo) {
      var v = cellInfo.value, m = cellInfo.row['%s']
      if (v === null || v === undefined || v === '') return ''
      var stat = cellInfo.row['depth_stat'] ? cellInfo.row['depth_stat'] : 'median'
      var mean = cellInfo.row['depth_mean']
      var tip = stat + ' depth ' + v + 'x'
      if (mean !== null && mean !== undefined && mean !== '' && stat === 'median') tip += ' (mean ' + mean + 'x)'
      if (!m || m <= 0) return '<span title=\"' + tip + '\">' + v + '</span>'
      var f = Math.max(0, Math.min(1, v / m))
      var w = Math.round(f * 46)
      var col = f < 0.05 ? '#c0392b' : (f < 0.33 ? '#b8860b' : '#2d6a4f')
      tip += ', ' + Math.round(f*100) + '%% of this sample&#39;s deepest scaffold'
      if (v === 0) tip += ' -- over half this scaffold has no coverage'
      return '<span title=\"' + tip + '\">' +
        '<span style=\"display:inline-block; width:48px; height:8px; vertical-align:middle;\">' +
        '<span style=\"display:inline-block; width:' + w + 'px; height:8px; background:' + col + ';\"></span></span>' +
        '<span style=\"margin-left:5px; font-size:0.85em; color:#555;\">' + v + '</span></span>'
    }", max_col) |> htmlwidgets::JS()
}

#' Render the five rDNA elements as filled/hollow pills
#' @noRd
rt_element_strip <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var want = ['18S','ITS1','5.8S','ITS2','28S']
      var have = (cellInfo.value ? String(cellInfo.value) : '').split(',')
      var out = ''
      for (var i = 0; i < want.length; i++) {
        var on = have.indexOf(want[i]) >= 0
        out += '<span style=\"display:inline-block; margin-right:2px; padding:0 4px; font-size:0.72em;' +
          ' border-radius:2px; border:1px solid ' + (on ? '#2d6a4f' : '#ccc') + '; color:' +
          (on ? '#fff' : '#bbb') + '; background:' + (on ? '#2d6a4f' : 'transparent') + ';\">' +
          want[i] + '</span>'
      }
      return '<span title=\"rRNA/ITS features annotated on this scaffold\">' + out + '</span>'
    }")
}

#' Format a unix timestamp as a date
#' @noRd
rt_ts_date <- function() {
  htmlwidgets::JS("
    function(cellInfo) {
      var options = { year: 'numeric', month: 'numeric', day: 'numeric', hour: 'numeric', minute: 'numeric' };
      var date = new Date(1000*cellInfo.value).toLocaleDateString('en-US', options);
      return date !== 'Invalid Date' ? date : null;
    }")
}

#' Cell: taxon name as a hyperlink to its NCBI Datasets taxonomy page.
#' Reads the row's `taxid`; shows "unassigned" (red) when there is no name.
#' @noRd
rt_taxid_link <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var name = cellInfo.value;
      var taxid = cellInfo.row ? cellInfo.row['taxid'] : null;
      if (!name) return '<span style=\"color:#c0392b;\">unassigned</span>';
      if (!taxid) return name;
      var url = 'https://www.ncbi.nlm.nih.gov/datasets/taxonomy/' + taxid + '/';
      return `<a href='${url}' target='_blank' rel='noopener'>${name} <span style='color:#888;'>(${taxid})</span></a>`;
    }")
}

#' Cell: taxonomy match confidence as a colored badge.
#' @noRd
rt_confidence_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (!v) v = 'unresolved';
      var colors = {exact:'#2d6a4f', manual:'#1f6feb', fuzzy_species:'#b8860b',
                    fuzzy_genus:'#b8860b', unresolved:'#c0392b'};
      var c = colors[v] || '#888';
      return `<span style='background:${c}; color:#fff; padding:1px 7px; border-radius:10px; font-size:0.8em; white-space:nowrap;'>${v}</span>`;
    }")
}

#' Cell: a clickable Font Awesome icon that reports the row (like rt_link but an
#' icon instead of the cell value). Used for the Assemble Details/Output buttons.
#' @noRd
rt_icon_link <- function(InputId, icon, title = "") {
  sprintf(
    "function(cellInfo) {
      var clickid = '%s';
      var sampid = cellInfo.index + 1;
      return `<a href='#' id=${sampid} class='grow' title='%s' ` +
        `onclick='event.stopPropagation(); Shiny.setInputValue(&#39;${clickid}&#39;, this.id, {priority: &#39;event&#39;})'>` +
        `<i class='%s'></i></a>`;
    }", InputId, title, icon
  ) |> htmlwidgets::JS()
}

#' Inline depth sparkline from a comma-separated track
#'
#' Shown because a curator deciding keep/drop should see coverage shape, and
#' deliberately not used by any rule: high depth is often k-mer collapse across
#' the conserved core or IGS subrepeats, and low depth is often isomer paths
#' splitting reads.
#' @noRd
rt_depth_sparkline <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (v === null || v === undefined || v === '' || v === 'NA') return '';
      var a = String(v).split(',').map(Number).filter(function(x){ return !isNaN(x); });
      if (a.length < 2) return '';
      var w = 90, h = 18, mx = Math.max.apply(null, a);
      if (!(mx > 0)) return '';
      var pts = a.map(function(y, i) {
        var px = (i / (a.length - 1)) * w;
        var py = h - (y / mx) * h;
        return px.toFixed(1) + ',' + py.toFixed(1);
      }).join(' ');
      return '<svg width=\"' + w + '\" height=\"' + h + '\" ' +
        'viewBox=\"0 0 ' + w + ' ' + h + '\" preserveAspectRatio=\"none\">' +
        '<polyline points=\"' + pts + '\" fill=\"none\" stroke=\"#6c757d\" ' +
        'stroke-width=\"1\"/></svg>';
    }")
}

#' Row click -> "examine this unit", for tables whose tick boxes mean something else
#'
#' The scaffold table's tick boxes feed Copy FASTA and bulk K/D, so the row click is free to
#' mean "show this scaffold in the panels below". Clicks on the tick box and on
#' the in-cell action buttons are ignored: the buttons already stop propagation,
#' and `.selection` would otherwise both tick the row and move the panels.
#' @param InputId Shiny input to set to the clicked row's uid.
#' @noRd
rt_examine_row <- function(InputId) {
  sprintf(
    "function(rowInfo, column) {
      if (column.id === '.selection' || column.id === 'triage') return;
      var uid = rowInfo.row['uid'];
      if (uid === null || uid === undefined) return;
      Shiny.setInputValue('%s', uid, {priority: 'event'});
    }", InputId) |> htmlwidgets::JS()
}

#' Row style for the scaffold table: grey when dropped, marked when examined
#'
#' The examined unit is read from the table's `meta`, which
#' [reactable::updateReactable()] can change without a re-render, so clicking a
#' row does not cost the user their sort, filters and page.
#' @noRd
rt_scaffold_row_style <- function() {
  htmlwidgets::JS("
    function(rowInfo, state) {
      var st = {};
      if (rowInfo.row['ignore'] === 1) { st.background = '#f3f3f3'; st.color = '#999'; }
      if (state.meta && state.meta.examined === rowInfo.row['uid']) {
        st.background = '#cadff5';
        st.boxShadow = 'inset 3px 0 0 #3c6fb5';
      }
      return st;
    }")
}
