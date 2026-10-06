window.MP_SAMPLE_TABLES = '#tbl';

// The main table's rows in display order, after its search and column filters.
window.rpVisibleRows = function() {
  try { return Reactable.getState('tbl').sortedData; } catch (e) { return []; }
};
// Report those rows (1-based data indices) whenever the table redraws. Only
// #tbl is watched, and reporting changes nothing in it, so an idle page does
// no work.
$( document ).ready(function(){
  var tbl = document.getElementById('tbl'), pending = null;
  if (!tbl) return;
  function report() {
    try {
      Shiny.setInputValue('tbl_visible',
        {rows: Reactable.getInstance('tbl').rows.map(function(r) { return r.index + 1; })});
    } catch (e) {}
  }
  function queue() { clearTimeout(pending); pending = setTimeout(report, 50); }
  new MutationObserver(queue).observe(tbl, { childList: true, subtree: true, characterData: true });
  // a data update can move rows off the page without changing the page itself
  Shiny.addCustomMessageHandler('rp_tbl_report', queue);
});

// Copy an element's text (the run report's Copy button). Plain http off
// localhost has no Clipboard API, so select the text and copy that instead.
window.rpCopyEl = function(id) {
  var el = document.getElementById(id), t = el.innerText;
  if (navigator.clipboard && window.isSecureContext) { navigator.clipboard.writeText(t); return; }
  var r = document.createRange(); r.selectNodeContents(el);
  var s = window.getSelection(); s.removeAllRanges(); s.addRange(r);
  document.execCommand('copy'); s.removeAllRanges();
};


// Sample-table resize: the CSS handle (resize: vertical) writes an inline
// height, which then sticks and stops the table re-fitting the window. On
// window resize, clear that inline height so the stylesheet calc() takes
// over again and the table stays fully displayed.
$( document ).ready(function(){
  window.addEventListener('resize', function() {
    document.querySelectorAll('.mp-table-resize').forEach(function(el) {
      el.style.height = '';
    });
  });
});

// The alignment viewer uses native scrolling only: a vertical wheel scrolls
// the page, horizontal input scrolls the view.

// msa.js measures its container when it renders. A modal that is still fading
// in can measure as zero width, and msa.js then falls back to the page width
// and overflows the dialog. Once the modal is shown, ask the server to draw
// again any alignment that overflowed. Same when a collapsed section opens.
$( document ).ready(function(){
  function relayout(root) {
    var over = Array.prototype.some.call(root.querySelectorAll('.msaR'), function(el) {
      return el.parentNode && el.scrollWidth > el.parentNode.clientWidth + 2;
    });
    if (over) Shiny.setInputValue('msa_relayout', Date.now(), {priority: 'event'});
  }
  $(document).on('shown.bs.modal', '#shiny-modal', function() { relayout(this); });
  $(document).on('hidden.bs.modal', '#shiny-modal', function() {
    Shiny.setInputValue('modal_hidden', Date.now(), {priority: 'event'});
  });
  document.addEventListener('toggle', function(e) {
    var d = e.target;
    if (d.tagName !== 'DETAILS' || !d.open) return;
    $(d).trigger('shown');
    setTimeout(function() { relayout(d); }, 0);
  }, true);
});

// Add an "All" choice to the sample-table page-size dropdowns. reactable
// (0.4.5) has no native "All", so append an option with a very large page size
// (shows every filtered row). A MutationObserver re-adds it after reactable
// re-renders the select.
$( document ).ready(function(){
  var GATED_ALL = window.MP_SAMPLE_TABLES;
  var ALL_PAGE_SIZE = 1000000;

  function addAllOption(select) {
    if (!select || select.querySelector('option[data-mp-all]')) return;
    var opt = document.createElement('option');
    opt.value = String(ALL_PAGE_SIZE);
    opt.text = 'All';
    opt.setAttribute('data-mp-all', '1');
    select.appendChild(opt);
  }

  function refreshAllOptions() {
    document.querySelectorAll(GATED_ALL).forEach(function(tbl) {
      tbl.querySelectorAll('.rt-page-size-select').forEach(addAllOption);
    });
  }

  refreshAllOptions();
  new MutationObserver(refreshAllOptions)
    .observe(document.body, { childList: true, subtree: true });
});

// Shift-click range selection for the main sample reactable tables.
// reactable (0.4.5) has no native range selection, so we drive it on the
// client: a plain click stores an "anchor" row, and a shift-click selects
// every row between the anchor and the clicked row by programmatically
// clicking their selection checkboxes. Additive only (never deselects).
$( document ).ready(function(){
  var GATED = window.MP_SAMPLE_TABLES;
  var anchors = {}; // per-table anchor index, keyed by table element id

  // Rows on the current page that carry a selection checkbox, in visual
  // (DOM) order. Excludes the select-all header (in .rt-thead) and any
  // pagination padding rows (no checkbox).
  function selectableRows(tbody) {
    return Array.prototype.filter.call(
      tbody.querySelectorAll('.rt-tr'),
      function(r) { return r.querySelector('input.rt-select-input'); }
    );
  }

  // Stop shift-click from highlighting cell text inside the gated tables.
  document.addEventListener('mousedown', function(e) {
    if (e.shiftKey && e.target.closest && e.target.closest(GATED)) {
      e.preventDefault();
    }
  }, true);

  document.addEventListener('click', function(e) {
    if (!e.isTrusted) return;                 // ignore our own synthetic clicks
    if (!e.target.closest) return;
    if (e.target.closest('a, button')) return; // leave link/icon actions alone
    var table = e.target.closest(GATED);
    if (!table) return;
    var row = e.target.closest('.rt-tr');
    if (!row || !row.querySelector('input.rt-select-input')) return;
    var tbody = row.closest('.rt-tbody');
    if (!tbody) return;

    var rows = selectableRows(tbody);
    var idx = rows.indexOf(row);
    if (idx < 0) return;
    var id = table.id;
    var anchor = anchors[id];

    if (e.shiftKey && anchor != null && anchor < rows.length) {
      var lo = Math.min(anchor, idx);
      var hi = Math.max(anchor, idx);
      for (var i = lo; i <= hi; i++) {
        if (i === idx) continue;              // reactable toggles the clicked row
        var cb = rows[i].querySelector('input.rt-select-input');
        if (cb && !cb.checked) cb.click();    // additive: only add unchecked rows
      }
      if (window.getSelection) window.getSelection().removeAllRanges();
      // Keep the anchor so successive shift-clicks re-extend from it.
    } else {
      anchors[id] = idx;                      // plain click sets a new anchor
    }
  }, true);
});

// Name the current step in the browser tab.
$(document).on('shiny:inputchanged', function(e) {
  if (e.name === 'mode' && e.value) document.title = 'RiboPilot - ' + e.value;
});
$(document).on('shiny:connected', function() {
  document.title = 'RiboPilot - ' + ($('#mode input:checked').val() || $('#mode').val() || 'Taxonomy');
});

// Export Data token chips: a chip inserts its text at the cursor of the header
// box last focused in the same modal (the FASTA header box by default).
$(document).on('focusin', 'textarea', function() {
  var list = $(this).closest('.modal').find('.mp-token-list');
  if (list.length) list.attr('data-target', this.id);
});
$(document).on('click', '.mp-token-chip', function(e) {
  e.preventDefault();
  var box = document.getElementById($(this).closest('.mp-token-list').attr('data-target'));
  if (!box) return;
  var ins = this.getAttribute('data-insert');
  var s = box.selectionStart, t = box.selectionEnd, v = box.value;
  box.value = v.slice(0, s) + ins + v.slice(t);
  box.selectionStart = box.selectionEnd = s + ins.length;
  box.focus();
  $(box).trigger('input').trigger('change');
});
// Chips whose column is empty for some records of the chosen export group turn
// warning-orange, with the count in the tooltip.
function mpTokenFlags(modal) {
  var list = modal.find('.mp-token-list');
  if (!list.length) return;
  var group = modal.find('select[id$="export_group_sel"]').val();
  var totals = JSON.parse(list.attr('data-totals') || '{}');
  var n = totals[group] || 0;
  list.find('.mp-token-chip').each(function() {
    var m = JSON.parse(this.getAttribute('data-missing') || '{}')[group] || 0;
    var base = this.getAttribute('data-title') || '';
    $(this).toggleClass('mp-token-missing', m > 0);
    this.title = m > 0 ? base + '\nMissing for ' + m + ' of ' + n + ' records in this group' : base;
  });
}
$(document).on('change', 'select[id$="export_group_sel"]', function() {
  mpTokenFlags($(this).closest('.modal'));
});
$(document).on('shown.bs.modal', function(e) { mpTokenFlags($(e.target).closest('.modal')); });
$(document).on('input', '.mp-token-filter', function() {
  var q = this.value.toLowerCase();
  var list = $(this).closest('.mp-token-list');
  list.find('.mp-token-chip').each(function() {
    $(this).toggle($(this).text().toLowerCase().indexOf(q) !== -1);
  });
  if (q) list.find('details').attr('open', '');
});

// "Insert field" beside a header box: aim the chips at that box, open the panel
$(document).on('click', '.mp-insert-field', function(e) {
  e.preventDefault();
  var modal = $(this).closest('.modal');
  var panel = modal.find('details.mp-meta-panel');
  if (!panel.length) return;
  modal.find('.mp-token-list').attr('data-target', this.getAttribute('data-box'));
  panel.prop('open', true);
  panel[0].scrollIntoView({behavior: 'smooth', block: 'start'});
  var f = panel.find('.mp-token-filter')[0];
  if (f) f.focus({preventScroll: true});
});
