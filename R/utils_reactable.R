# reactable rendering helpers (mirrors MitoPilot's utils_reactable.R).

#' Highlight selected row(s)
#' @noRd
rt_highlight_row <- function() {
  htmlwidgets::JS("
    function(rowInfo) {
      if (typeof rowInfo === 'undefined') return
      var col = rowInfo.selected ? '#D3BEC2' : '#FFFFFF'
      return { background: col }
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
