# reactable rendering helpers (mirrors MitoPilot's utils_reactable.R).

#' Highlight selected row(s)
#'
#' @noRd
rt_highlight_row <- function() {
  htmlwidgets::JS(
    "
    function(rowInfo) {
      if ( typeof rowInfo === 'undefined') return
      var col = rowInfo.selected ? 'var(--mp-primary-soft)' : 'var(--mp-surface)'
      return { background: col }
    }
    "
  )
}

#' Escape a label for interpolation into a single-quoted HTML attribute
#'
#' Every helper below builds HTML inside a JS template literal, so one
#' apostrophe or angle bracket in a caller-supplied label silently breaks the
#' cell.
#'
#' @noRd
mp_js_attr <- function(x) {
  htmltools::htmlEscape(as.character(x), attribute = TRUE)
}

#' Build a JS object literal from a named vector, values escaped
#'
#' @noRd
mp_js_obj <- function(x) {
  if (length(x) == 0) {
    return("{}")
  }
  paste0(
    "{ ",
    paste(sprintf("'%s': '%s'", names(x), mp_js_attr(x)), collapse = ", "),
    " }"
  )
}

#' Lock glyph shared by the main table lock columns
#' @noRd
rt_lock_icon <- function(inputId = NULL) {
  rt_dynamicIcon(c(`0` = "fa fa-lock-open mp-fg-neutral", `1` = "fa fa-lock mp-fg-neutral"),
    labels = c(`0` = "Unlocked", `1` = "Locked"), inputId = inputId)
}

#' Dynamic Icon
#'
#' @param icons a named character vector list of icons to use.
#' @param labels optional named character vector, keyed identically to
#'   `icons`, giving the accessible name of each state. When supplied the
#'   glyph gets `role="img"`, `aria-label` and `title`; without it the markup
#'   is unchanged.
#' @param inputId optional Shiny input id. When supplied the icon is a button
#'   that sends its 1-based row index to this input.
#'
#' @noRd
rt_dynamicIcon <- function(icons = NULL, labels = NULL, inputId = NULL) {
  if (length(icons) == 0) {
    return({
      htmlwidgets::JS("function(cellInfo) {return cellInfo.value}")
    })
  }

  sprintf(
    "
    function(cellInfo){
      var { value } = cellInfo;
      var icons = %s;
      var labels = %s;
      var icon = icons[value] || '';
      var label = labels[value];
      var name = label ?
        ` role='img' aria-label='${label}' title='${label}'` : '';
      var html = `<i class='${icon}'${name} ` +
        `style='padding-left: 0.2em;'></i>`;
      var clickid = '%s';
      if (!clickid) return html;
      return `<a href='#' class='grow mp-icon-btn' data-row='${cellInfo.index + 1}' ` +
        `onclick='event.preventDefault(); event.stopPropagation(); ` +
        `Shiny.setInputValue(&#39;${clickid}&#39;, this.dataset.row, {priority: &#39;event&#39;})'>` +
        html + `</a>`;
    }
    ",
    mp_js_obj(icons), mp_js_obj(labels), inputId %||% ""
  ) |> htmlwidgets::JS()
}

#' Column header with a tooltip (wrapper kept for existing call sites)
#' @noRd
rt_help_header <- function(name, tip) as.character(rt_header(name, tip))

#' Add hover info to truncated text in reactable
#'
#' Assumes that wrap=FALSE in the table so that long text is truncated.
#'
#' @noRd
rt_longtext <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var raw = (cellInfo.value === null || cellInfo.value === undefined) ?
        '' : String(cellInfo.value)
      if (raw === '') return ''
      var text = raw.replace(/&/g, '&amp;').replace(/</g, '&lt;')
        .replace(/'/g, '&#39;')
      return `<abbr class='mp-abbr' tabindex='0' ` +
      `title='${text}'>${text}</abbr>`
    }"
  )
}


#' Render a cell value as a status pill
#'
#' The one pill component. Colour carries the judgement, so the caller maps
#' each value to a meaning rather than to a colour.
#'
#' @param map named character vector, cell value -> one of "success",
#'   "warning", "danger", "neutral", "info".
#' @param labels optional named character vector, cell value -> display text.
#'   Defaults to the value itself.
#' @param empty text for a missing/blank value ("not set"). Use "" to render
#'   nothing at all.
#' @param hide character vector of values that render nothing.
#' @noRd
rt_pill <- function(map, labels = NULL, empty = "not set", hide = NULL) {
  hide_js <- if (length(hide) == 0) {
    "[]"
  } else {
    paste0("[", paste(sprintf("'%s'", mp_js_attr(hide)), collapse = ", "), "]")
  }
  sprintf(
    "function(cellInfo) {
      var map = %s;
      var labels = %s;
      var hide = %s;
      var empty = '%s';
      var v = cellInfo.value;
      var val = (v === null || v === undefined) ? '' : String(v);
      if (hide.indexOf(val) >= 0) return '';
      var tone, text;
      if (val === '') {
        if (empty === '') return '';
        tone = 'neutral';
        text = empty;
      } else {
        tone = map[val] || 'neutral';
        text = labels[val] || val;
      }
      return `<span class='mp-pill mp-pill-${tone}'>${text}</span>`;
    }",
    mp_js_obj(map), mp_js_obj(labels), hide_js, mp_js_attr(empty)
  ) |> htmlwidgets::JS()
}

#' Render a yes/no text column as a status pill
#'
#' "no" is grey in every column, never green; a blank value reads "not set".
#'
#' @param invert if TRUE, "yes" is the bad outcome and gets the red pill.
#' @param hide_no if TRUE, render an empty cell for "no"/NA values (only
#'   "yes" gets a pill). Useful for columns where "no" is the default
#'   and noisy to display.
#' @param yes_tone pill tone for "yes" when invert is FALSE; "warning" for a
#'   flag whose "yes" asks for attention (problematic, partial).
#' @noRd
rt_bool_badge <- function(invert = FALSE, hide_no = FALSE, yes_tone = "success", no_tone = "neutral") {
  rt_pill(
    map = c(yes = if (invert) "danger" else yes_tone, no = no_tone),
    empty = if (hide_no) "" else "not set",
    hide = if (hide_no) "no" else NULL
  )
}

#' Render a curation flag count: blank/green when 0, orange badge when > 0
#' @noRd
rt_flag_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var n = cellInfo.value
      if (n === null || n === undefined || n === '' ) return ''
      n = Number(n)
      if (n === 0) return '<span class=\"mp-fg-success\">&#10003;</span>'
      var tip = cellInfo.row.curate_flags || ''
      return '<span title=\"' + String(tip).replace(/\"/g, '&quot;') + '\" class=\"mp-pill mp-pill-warning\">' + n + '</span>'
    }")
}

#' Render the stock-ITSx check: tick when clean, grey when not run, amber pill
#' with the flag list otherwise (or the 'n of m units' roll-up on a sample row)
#' @noRd
rt_itsx_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === null || v === undefined || v === '') return ''
      v = String(v)
      var esc = v.replace(/\"/g, '&quot;')
      if (v.indexOf('not_run(') === 0) return '<span title=\"' + esc + '\" style=\"color:#999;\">not run</span>'
      if (/^0 of /.test(v)) return '<span class=\"mp-fg-success\">&#10003;</span>'
      var txt = v.split(';').join('; ')
      return '<span title=\"' + esc + '\" class=\"mp-pill mp-pill-warning\">' + txt + '</span>'
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
      if (n === 0) return '<span class=\"mp-pill mp-pill-danger\">0</span>'
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

#' Add text click action to a cell
#'
#' @param InputId shiny input id to use
#' @param title optional tooltip naming the action. Without it the tooltip
#'   echoes the cell value, which is noise.
#' @noRd
rt_link <- function(InputId, title = NULL, lock_col = NULL) {
  tip <- if (is.null(title)) "${cellInfo.value}" else mp_js_attr(title)
  # On a locked row the cell is plain text: the click falls through to row
  # selection instead of silently doing nothing (theme T01).
  locked <- if (is.null(lock_col)) "" else sprintf(
    "var lk = (cellInfo.row || {})['%s'];
                if (lk == 1) {
                  return `<span class='mp-locked-cell' title='Locked - options are set to ${cellInfo.value}. Select the row and press Unlock to change them.'>${cellInfo.value}</span>`;
                }",
    lock_col
  )
  sprintf(
    "function(cellInfo) {
                // An empty cell is not a link: nothing to click through to.
                if (cellInfo.value === null || cellInfo.value === undefined ||
                    cellInfo.value === '') { return ''; }
                %s
                var clickid = '%s';
                var sampid = cellInfo.index+1;
                return `<a href='#' id='${clickid}_${sampid}' class='grow' title='%s' ` +
                `onclick='event.stopPropagation(); Shiny.setInputValue(&#39;${clickid}&#39;, ${sampid}, {priority: &#39;event&#39;})'>` +
                cellInfo.value +
                `</a>`;
                }",
    locked, InputId, tip
  ) |>
    htmlwidgets::JS()
}

#' Column header with a tooltip
#'
#' Pass to `colDef(header = )` alongside the usual `name`, so the header keeps
#' its accessible name and gains a one-sentence explanation on hover.
#'
#' @param name header text
#' @param tip one sentence explaining the column. NULL returns plain text.
#' @noRd
rt_header <- function(name, tip = NULL) {
  if (is.null(tip) || !nzchar(tip)) {
    return(name)
  }
  htmltools::tags$span(class = "mp-th-tip", title = tip, name)
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
      var tone
      if (close[v])      { tone = 'success' }
      else if (mid[v])   { tone = 'success' }
      else if (v === 'none') { tone = 'danger' }
      else               { tone = 'warning' }
      var t = cellInfo.row['screen_species'] ? cellInfo.row['screen_species'] : ''
      return '<span title=\"top BLAST hit: ' + t + '\" class=\"mp-pill mp-pill-' + tone + '\">' + v + '</span>'
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
        return '<span title=\"BLAST at NCBI core_nt\" class=\"mp-pill mp-pill-success\">NCBI</span>'
      if (v === 'local_rrna')
        return '<span title=\"BLAST against the bundled rDNA database\"' +
          ' class=\"mp-pill mp-pill-success\">bundled</span>'
      if (v === 'seeds')
        return '<span title=\"local BLAST against the GetOrganelle seed file: the closest bait' +
          ' sequence, NOT a GenBank identification\" class=\"mp-pill mp-pill-warning\">seed file</span>'
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
        sib + ' at ' + rk + ')\" class=\"mp-pill mp-pill-warning\">odd</span>'
    }")
}

#' Render the contamination-screen verdict: is the scaffold on this sample's
#' own lineage.
#'
#' `no_reference` renders NEUTRAL GREY like `rt_shares_badge()`'s "no
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
      var tone
      if (v === 'on_target')       { tone = 'success' }
      else if (v === 'off_target') { tone = 'danger' }
      else                          { tone = 'warning' }
      return '<span class=\"mp-pill mp-pill-' + tone + '\">' + v + '</span>'
    }")
}

#' Render the chimera check verdict: clean, minor, or severe, with the
#' reason as a hover tooltip.
#' @noRd
rt_chimera_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (v === null || v === undefined || v === '' || v === 'NA') return '';
      var tone;
      if (v === 'clean')        { tone = 'success'; }
      else if (v === 'severe')  { tone = 'danger'; }
      else                      { tone = 'warning'; }
      var tip = cellInfo.row.chimera_reason || '';
      return '<span title=\"' + String(tip).replace(/\"/g, '&quot;') + '\" class=\"mp-pill mp-pill-' + tone + '\">' + v + '</span>';
    }")
}

#' Render the Checks summary: the pill colour is the worst finding's tone
#' (`checks_tone`), and the text names every finding, so colour is never the
#' only signal
#' @noRd
rt_checks_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (v === null || v === undefined || v === '' || v === 'NA') return '';
      var tone = cellInfo.row.checks_tone || 'neutral';
      var t = String(v).replace(/\"/g, '&quot;');
      return '<span title=\"' + t + '\" class=\"mp-pill mp-pill-' + tone + '\"' +
        ' style=\"max-width:100%; overflow:hidden; text-overflow:ellipsis; vertical-align:middle;\">' + v + '</span>';
    }")
}

#' Render the intragenomic variant class: is this scaffold one of several rDNA
#' copies, or the same molecule assembled twice.
#'
#' NEITHER class is a defect, so neither may borrow `rt_screen_badge()`'s
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
#' The value is the MEDIAN per-position depth (see R/coverage.R for why a mean
#' misleads on an rDNA scaffold); the tooltip also shows the mean for context.
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
      var mean = cellInfo.row['depth_mean']
      var tip = 'median depth ' + v + 'x'
      if (mean !== null && mean !== undefined && mean !== '') tip += ' (mean ' + mean + 'x)'
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
      var date = new Date(1000*cellInfo.value).toLocaleDateString(undefined, options);
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
                    fuzzy_genus:'#b8860b', ambiguous:'#c0392b', unresolved:'#c0392b'};
      var c = colors[v] || '#888';
      return `<span style='background:${c}; color:#fff; padding:1px 7px; border-radius:10px; font-size:0.8em; white-space:nowrap;'>${v}</span>`;
    }")
}


#' Add reactable icon button with text
#'
#' @param inputId shiny input id to use
#' @param icon font awesome icon name
#' @param text fallback label used when the cell value is empty. The cell
#'   value is still the render gate: an empty cell renders no button.
#' @param label optional explicit button text. When supplied it replaces the
#'   cell value in the label only, so rows with nothing behind them still
#'   render no button.
#' @param title optional tooltip naming the action.
#'
#' @noRd
rt_icon_bttn_text <- function(inputId, icon, text = "", label = NULL,
                              title = NULL) {
  relabel <- if (is.null(label)) {
    ""
  } else {
    sprintf("value = '%s';", mp_js_attr(label))
  }
  tip <- if (is.null(title)) {
    ""
  } else {
    sprintf(" title='%s'", mp_js_attr(title))
  }
  sprintf(
    "
    function(cellInfo) {
      var { index, value } = cellInfo;
      value = value ? value : '%s';
      if (value === undefined || value === null || value==='') {
        return;
      }
      %s
      return `<button type='button' ` +
        `class='icon-bttn-text grow' ` +
        `id='%s_${index+1}'%s ` +
        `onclick='event.stopPropagation(); Shiny.setInputValue(&#39;%s&#39;, ${index+1}, {priority: &#39;event&#39;})'>` +
        `<i class='%s' aria-hidden='true' ` +
        `style='margin-right: 4px;'></i>` +
        `<small>${value}</small>` +
        `</button>`
    }
    ",
    mp_js_attr(text), relabel, inputId, tip, inputId, icon
  ) |>
    htmlwidgets::JS()
}

#' Inline depth sparkline from a comma-separated track
#'
#' Shown because a curator deciding keep/drop should see coverage shape, and
#' deliberately not used by any rule: high depth is often k-mer collapse across
#' the conserved core or IGS subrepeats, and low depth is often a repeat inside
#' the scaffold, where reads map to two places and are left out.
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
