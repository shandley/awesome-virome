// This function initializes mermaid
function initializeMermaid() {
  // Check if mermaid is loaded
  if (typeof window.mermaid !== 'undefined') {
    try {
      // Configure and initialize
      window.mermaid.initialize({
        startOnLoad: true,
        theme: 'default',
        flowchart: {
          useMaxWidth: false,
          htmlLabels: true,
          curve: 'basis'
        },
        securityLevel: 'loose'
      });
    } catch (e) {
      console.error('Mermaid initialization error:', e);
    }
  } else {
    console.log('Mermaid not defined, loading script...');
    // Load mermaid if not already loaded
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/mermaid@9.4.3/dist/mermaid.min.js';
    script.async = true;
    script.onload = function() {
      console.log('Mermaid loaded, initializing...');
      // Initialize after loading
      if (typeof window.mermaid !== 'undefined') {
        window.mermaid.initialize({
          startOnLoad: true,
          theme: 'default',
          flowchart: {
            useMaxWidth: false,
            htmlLabels: true,
            curve: 'basis'
          },
          securityLevel: 'loose'
        });
      }
    };
    document.head.appendChild(script);
  }
}

// Try to initialize as soon as possible
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', initializeMermaid);
} else {
  // DOM already loaded, initialize now
  initializeMermaid();
}

// Add a fallback to make sure it runs
window.addEventListener('load', function() {
  console.log('Window loaded, initializing mermaid if needed...');
  // Re-initialize diagrams that might have been missed
  if (typeof window.mermaid !== 'undefined') {
    try {
      window.mermaid.init(undefined, document.querySelectorAll('.mermaid'));
    } catch (e) {
      console.error('Mermaid re-initialization error:', e);
    }
  }
});