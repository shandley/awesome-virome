document.addEventListener('DOMContentLoaded', function() {
  if (typeof mermaid !== 'undefined') {
    mermaid.initialize({
      startOnLoad: true,
      theme: 'default',
      flowchart: {
        useMaxWidth: false,
        htmlLabels: true,
        curve: 'basis'
      },
      securityLevel: 'loose'
    });
  } else {
    console.log('Mermaid not loaded yet, loading now...');
    var script = document.createElement('script');
    script.onload = function() {
      mermaid.initialize({
        startOnLoad: true,
        theme: 'default',
        flowchart: {
          useMaxWidth: false,
          htmlLabels: true,
          curve: 'basis'
        },
        securityLevel: 'loose'
      });
    };
    script.src = 'https://cdn.jsdelivr.net/npm/mermaid@9.4.3/dist/mermaid.min.js';
    document.head.appendChild(script);
  }
});