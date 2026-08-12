/* Chunked / resumable file upload for the CRISPRme Settings page.
 *
 * Large genome/VCF files cannot go through Dash's dcc.Upload (base64, held in
 * memory). This slices the chosen file into fixed-size chunks and POSTs them to
 * the Flask endpoint /settings/upload-chunk one at a time, so browser memory use
 * stays flat regardless of file size. Widgets are plain <div> elements with
 * class "crisprme-chunk-upload" and data-target (+ optional data-dataset-input);
 * a MutationObserver wires them up as Dash renders them into the SPA.
 */
(function () {
  var CHUNK = 8 * 1024 * 1024; // 8 MB per chunk

  function uploadFile(file, target, dataset, genome, prog, btn) {
    var total = Math.max(1, Math.ceil(file.size / CHUNK));
    var idx = 0;
    btn.disabled = true;

    function sendNext() {
      if (idx >= total) return;
      var start = idx * CHUNK;
      var blob = file.slice(start, Math.min(start + CHUNK, file.size));
      fetch("/settings/upload-chunk", {
        method: "POST",
        headers: {
          "X-Target": target,
          "X-File-Name": file.name,
          "X-Dataset": dataset,
          "X-Genome": genome,
          "X-Chunk-Index": String(idx),
          "X-Total-Chunks": String(total)
        },
        body: blob
      })
        .then(function (r) {
          return r.text().then(function (t) {
            return { ok: r.ok, text: t };
          });
        })
        .then(function (res) {
          if (!res.ok) {
            prog.textContent = "Error: " + res.text;
            btn.disabled = false;
            return;
          }
          idx += 1;
          var pct = Math.round((idx / total) * 100);
          prog.textContent = "Uploading... " + pct + "%";
          if (res.text.indexOf("complete") === 0) {
            prog.textContent = "Upload complete. It appears in the search form and in Settings when you revisit the page.";
            btn.disabled = false;
            return;
          }
          sendNext();
        })
        .catch(function (e) {
          prog.textContent = "Upload failed: " + e;
          btn.disabled = false;
        });
    }
    sendNext();
  }

  function wire(widget) {
    if (widget.dataset.wired) return;
    widget.dataset.wired = "1";
    var target = widget.dataset.target || "";
    var datasetInputId = widget.dataset.datasetInput || "";
    var genomeInputId = widget.dataset.genomeInput || "";

    var input = document.createElement("input");
    input.type = "file";
    var btn = document.createElement("button");
    btn.type = "button";
    btn.textContent = "Upload file";
    btn.style.marginLeft = "8px";
    var prog = document.createElement("span");
    prog.style.marginLeft = "8px";
    prog.style.color = "#555";

    widget.appendChild(input);
    widget.appendChild(btn);
    widget.appendChild(prog);

    btn.addEventListener("click", function () {
      var file = input.files && input.files[0];
      if (!file) {
        prog.textContent = "Choose a file first.";
        return;
      }
      var dataset = "";
      if (datasetInputId) {
        var di = document.getElementById(datasetInputId);
        if (di) dataset = di.value || "";
      }
      var genome = "";
      if (genomeInputId) {
        // dcc.Dropdown renders its value in a hidden input inside the container
        var gc = document.getElementById(genomeInputId);
        if (gc) {
          var gi = gc.querySelector("input");
          if (gi && gi.value) genome = gi.value;
          var sv = gc.querySelector("[class*=singleValue]");
          if (!genome && sv) genome = sv.textContent || "";
        }
      }
      uploadFile(file, target, dataset, genome, prog, btn);
    });
  }

  function scan() {
    var widgets = document.querySelectorAll(".crisprme-chunk-upload");
    for (var i = 0; i < widgets.length; i++) wire(widgets[i]);
  }

  document.addEventListener("DOMContentLoaded", scan);
  var mo = new MutationObserver(scan);
  mo.observe(document.documentElement, { childList: true, subtree: true });
})();
