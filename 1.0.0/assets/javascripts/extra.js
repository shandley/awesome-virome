document.addEventListener("DOMContentLoaded", function() {
  // Wait until mermaid is fully loaded
  setTimeout(function() {
    if (typeof mermaid !== 'undefined') {
      mermaid.initialize({
        startOnLoad: true,
        theme: "default",
        flowchart: {
          useMaxWidth: false,
          htmlLabels: true,
          curve: 'basis'
        },
        securityLevel: 'loose'
      });
      // Force rendering of any mermaid diagrams
      mermaid.init(undefined, document.querySelectorAll('.mermaid'));
    }
  }, 500);
});