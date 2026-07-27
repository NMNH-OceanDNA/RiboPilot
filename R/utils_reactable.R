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
  key <- purrr::imap_chr(icons, ~ paste0(.y, ": '", .x, "'")) |> unname() |> paste(collapse = ", ")
  sprintf(
    "function(cellInfo){
       var { value } = cellInfo;
       var icons = { %s };
       var icon = icons[value] || '';
       return `<i class='${icon}' style='padding-left: 0.2em;'></i>`
     }", key
  ) |> htmlwidgets::JS()
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

#' Clickable cell keyed on the row's own `uid`, not its position
#'
#' [rt_link()] sends `cellInfo.index + 1`, so every handler has to re-query and
#' index back into the result: if the two orderings ever disagree the click acts
#' on the WRONG row, and on this table that means deleting another scaffold's
#' curation. The uid (`ID|path|scaffold`) is carried in the row itself, so
#' sorting, filtering and paging cannot desynchronise it. IDs are constrained to
#' `^[a-zA-Z0-9_:-]+$` at project creation, so `|` cannot occur inside one.
#' @param InputId shiny input id to set on click.
#' @param suffix Optional string appended after a second `|` (e.g. an action).
#' @noRd
rt_unit_link <- function(InputId, suffix = NULL) {
  sprintf(
    "function(cellInfo) {
      var clickid = '%s';
      var uid = cellInfo.row['uid'];
      if (uid === null || uid === undefined) return cellInfo.value ? cellInfo.value : '';
      var payload = uid%s;
      return `<a href='#' class='grow' ` +
        `onclick='event.stopPropagation(); Shiny.setInputValue(&#39;${clickid}&#39;, ` +
        `&#39;${payload}&#39;, {priority: &#39;event&#39;})'>` +
        (cellInfo.value ? cellInfo.value : '') + `</a>`;
    }", InputId, if (is.null(suffix)) "" else sprintf(" + '|%s'", suffix)
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
      function btn(label, on, colOn) {
        var st = on
          ? 'background:' + colOn + '; color:#fff; border:1px solid ' + colOn + ';'
          : 'background:transparent; color:#999; border:1px solid ' + (v ? '#ccc' : '#e0a800') + ';';
        return '<a href=\"#\" style=\"display:inline-block; width:20px; text-align:center; ' +
          'font-size:0.75em; line-height:16px; border-radius:3px; margin-right:2px; ' +
          'text-decoration:none; ' + st + '\" ' +
          'onclick=\"event.stopPropagation(); Shiny.setInputValue(&#39;' + clickid +
          '&#39;, &#39;' + uid + '|' + label.toLowerCase() + '&#39;, {priority: &#39;event&#39;})\">' +
          label + '</a>';
      }
      return '<span title=\"' + tip.replace(/\"/g, '') + '\">' +
        btn('K', keep, '#2d6a4f') + btn('D', drop, '#c0392b') + '</span>';
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

#' Inline depth bar relative to the deepest scaffold of the same sample
#'
#' The value is the MEDIAN per-position depth where one is stored (see
#' R/coverage.R for why a mean misleads on an rDNA scaffold); the tooltip names
#' which statistic it is, so a pre-v10 row falling back to the mean is not read
#' as a median.
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
        '<span style=\"display:inline-block; width:48px; height:8px; background:#eee; vertical-align:middle;\">' +
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

#' Cell: a greyed, non-clickable Font Awesome icon (a disabled button).
#' @noRd
rt_icon_disabled <- function(icon, title = "") {
  sprintf(
    "function(cellInfo) {
      return `<span title='%s' style='color:#bbb; cursor:not-allowed;'><i class='%s'></i></span>`;
    }", title, icon
  ) |> htmlwidgets::JS()
}
